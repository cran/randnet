#library(poweRlaw)
#library(RSpectra)
library(Matrix)
#library(irlba)
#library(AUC)
library(entropy)



BlockModel.Gen <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0,simple=TRUE,power=TRUE,alpha=5,degree.seed=NULL){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(rho*0.2+(1-rho))^2)
    if((rho >0) && (!simple) && (!power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(0.6)^2)   }
    if((rho >0) && (!simple) && (power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*((1.285)^2))   }

    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    for(i in 1:n){
        M[i,membership[i]] <- 1
    }
    A.bar <- M%*%P%*%t(M)
    node.degree <- rep(1,n)
    if(rho>0){
    if(simple){
        node.degree[runif(n)<rho] <- 0.2
    }else{
        if(power==FALSE){
            node.degree <- runif(n)*0.8 + 0.2
        }else{
            MM <- ceiling(n/300)
            if(is.null(degree.seed)){
            degree.seed <- rplcon(300,1,alpha)
            }### sample fixed number of values from power law and then randomly assign to be degrees. Easier to control noises across different sample sizes.
            node.degree <- sample(degree.seed,size=n,replace=TRUE)
        }
    }}
    DD <- diag(node.degree)
    A.bar <- DD%*%A.bar%*%DD
    A.bar <- A.bar*lambda/mean(colSums(A.bar))
    #diag(A.bar) <- 0
    #avg.d <- mean(colSums(A.bar))
    #A.bar <- A.bar*lambda/avg.d
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    diag(A) <- 0
    return(list(A=A,g=membership,P=A.bar,theta=node.degree))
}



reg.SP <- function(A,K,tau=1,lap=FALSE){
    avg.d <- mean(colSums(A))
    A.tau <- A + tau*avg.d/nrow(A)
    if(!lap){SVD <- irlba(A.tau,nu=K,nv=K)}else{
         d.tau <- colSums(A.tau)
         L.tau <- diag(1/sqrt(d.tau))%*%A.tau%*%diag(1/sqrt(d.tau))
         #SVD <- svd(L.tau,nu=K,nv=K)
         SVD <- irlba(L.tau,nu=K,nv=K)
    }
    km <- kmeans(SVD$v[,1:K],centers=K,nstart=30,iter.max=30)#,algorithm="Lloyd")
    return(list(cluster=km$cluster,loss=km$tot.withinss))
}


reg.SSP <- function(A,K,tau=1,lap=FALSE){
    avg.d <- mean(colSums(A))
    A.tau <- A + tau*avg.d/nrow(A)
    if(!lap){SVD <- irlba(A.tau,nu=K,nv=K)
         V <- SVD$v[,1:K]
         V.norm <- apply(V,1,function(x)sqrt(sum(x^2)))
         V.normalized <- diag(1/V.norm)%*%V
         }else{
         d.tau <- colSums(A.tau)
         L.tau <- diag(1/sqrt(d.tau))%*%A.tau%*%diag(1/sqrt(d.tau))
         #SVD <- svd(L.tau,nu=K,nv=K)
         SVD <- irlba(L.tau,nu=K,nv=K)
         V <- SVD$v[,1:K]
         V.norm <- apply(V,1,function(x)sqrt(sum(x^2)))
         V.normalized <- diag(1/V.norm)%*%V
    }
    km <- kmeans(V.normalized,centers=K,nstart=50,iter.max=50)#,algorithm="Lloyd")
    return(list(cluster=km$cluster,loss=km$tot.withinss))
}


SBM.estimate <- function(A,g){
    n <- nrow(A)
    K <- length(unique(g))
    B <- matrix(0,K,K)
    M <- matrix(0,n,K)
    for(i in 1:K){
        for(j in i:K){
            if(i!=j){
                B[j,i] <- B[i,j] <- mean(A[which(g==i),which(g==j)])
            }else{
                n.i <- length(which(g==i))
                B[i,i] <- sum(A[which(g==i),which(g==i)])/(n.i^2 - n.i)
            }
        }
    }
    M[matrix(c(1:n,g),ncol=2)] <- 1
    P <- M%*%B%*%t(M)
    return(list(B=B,Phat=P,g=g))
}



DCSBM.estimate <- function(A,g){
    K <- length(unique(g))
    n <- nrow(A)
        B <- matrix(0,K,K)
        Psi <- rep(0,n)
        #Psi.outer <- outer(Psi,Psi)
        Theta <- matrix(0,n,K)
        for(i in 1:K){
            for(j in i:K){
                N.i <- which(g==i)
                N.j <- which(g==j)
                B[i,j] <- B[j,i] <- sum(A[N.i,N.j],na.rm=TRUE)+1e-3
            }
            Theta[N.i,i] <- 1
        }

        Psi <- colSums(A,na.rm=TRUE)
            #print(Psi)
        B.rowSums <- rowSums(B)
        B.g <- Theta%*%B.rowSums
        Psi <- as.numeric(Psi/B.g)
            #print(Psi)
        tmp.mat <- Theta*Psi
        P.hat <- tmp.mat%*%B%*%t(tmp.mat)
    return(list(Phat=P.hat,B=B,Psi=Psi,g=g))

}


NMI <- function(g1,g2){
    cross.table <- table(g1,g2)
    return(mi.empirical(as.matrix(cross.table))/entropy.empirical(as.vector(cross.table)))
}



### Bickel and Wang 2015 method
LRBIC <- function(A,Kmax,lambda=NULL,model="both"){
    n <- nrow(A)
    if(is.null(lambda)) lambda <- 1/n
    d <- colSums(A)
    SBM.result <- DCSBM.result <- rep(0,Kmax)
    sbm.p <- sum(A)/(n*(n-1))
    upper.index <- which(upper.tri(A))
    a <- A[upper.index]
    sbm.ll <- sum(a*log(sbm.p)) + sum((1-a)*log(1-sbm.p))
    dc.P <- outer(d,d)/sum(A)
    dc.p <- dc.P[upper.index]
    dc.p[dc.p>(1-1e-6)] <- 1-1e-6
    dc.p[dc.p<1e-6] <- 1e-6
    dc.ll <- sum(a*log(dc.p)) + sum((1-a)*log(1-dc.p))
    SBM.result[1] <- sbm.ll-lambda*n*log(n)
    DCSBM.result[1] <- dc.ll-lambda*n*log(n)
    try.model <- model
    for(K in 2:Kmax){
        if(try.model=="both"){
            model <- "SBM"
        }
        ##  SBM
        if(model=="SBM"){
        SBM.clust <- reg.SP(A,K=K,tau=0.25,lap=TRUE)
        #DCSBM.clust <- Arash.reg.SSP(A,K=K,tau=0.25,lap=TRUE)
        g <- SBM.clust$cluster
        n.K <- as.numeric(table(g))
        Pi <- n.K/n
        B <- matrix(0,K,K)
        for(j in 1:(K)){
            for(k in j:K){
                j.index <- which(g==j)
                k.index <- which(g==k)
                if(j!=k){
                    B[j,k] <- mean(A[j.index,k.index])
                    B[k,j] <- B[j,k]
                }else{
                    B[j,j] <- sum(A[j.index,j.index])/(length(j.index)^2-length(j.index))
                }
            }
        }
        Z <- matrix(0,n,K)
        Z[cbind(1:n,g)] <- 1
        SBM.P <- Z%*%B%*%t(Z)
        sbm.p <- SBM.P[upper.index]
        sbm.p[sbm.p>(1-1e-6)] <- 1-1e-6
        sbm.p[sbm.p<1e-6] <- 1e-6
        sbm.ll <- sum(a*log(sbm.p)) + sum((1-a)*log(1-sbm.p)) + sum(n.K*log(Pi))
        SBM.result[K] <- sbm.ll-lambda*(K)*(K+1)*n*log(n)/2
    }
        ## DCSBM
        if(try.model=="both"){
            model=="DCSBM"
        }
        if(model=="DCSBM"){
        DCSBM.clust <- reg.SSP(A,K=K,tau=0.25,lap=TRUE)
        g <- DCSBM.clust$cluster
        n.K <- as.numeric(table(g))
        Pi <- n.K/n
        B.star <- matrix(0,K,K)

        for(j in 1:(K)){
            for(k in j:K){
                j.index <- which(g==j)
                k.index <- which(g==k)
                B.star[j,k] <- sum(A[j.index,k.index])
                B.star[k,j] <- B.star[j,k]
            }
        }
        b.row <- rowSums(B.star)
        Z <- matrix(0,n,K)
        Z[cbind(1:n,g)] <- 1
        Z.b.row <- Z%*%b.row
        theta <- d/Z.b.row
        Z.theta <- Z*as.numeric(theta)
        dc.P <- Z.theta%*%B.star%*%t(Z.theta)
        dc.p <- dc.P[upper.index]
        dc.p[dc.p>(1-1e-6)] <- 1-1e-6
        dc.p[dc.p<1e-6] <- 1e-6
        dcsbm.ll <- sum(a*log(dc.p)) + sum((1-a)*log(1-dc.p)) + sum(n.K*log(Pi))
        DCSBM.result[K] <- dcsbm.ll-lambda*(K)*(K+1)*n*log(n)/2
        }
    }
    SBM.K <- DCSBM.K <- NA
    if((try.model=="SBM")||(try.model=="both")){
        SBM.K <- which.max(SBM.result)
    }
    if((model=="DCSBM")||(try.model=="both")){
    DCSBM.K <- which.max(DCSBM.result)
    }
    return(list(SBM.K=SBM.K,DCSBM.K=DCSBM.K,SBM.BIC=SBM.result,DCSBM.BIC=DCSBM.result))
}






#library(RSpectra)


## Le & Levina 2015 method
BHMC.estimate <- function(A,K.max=15){
    if(K.max <= 2) K.max <- 2
    d <- colSums(A)
    n <- nrow(A)
    I <- as(diag(rep(1,n)),"dgCMatrix")
    D <- as(diag(d),"dgCMatrix")
    r <- sqrt(mean(d))#sqrt(sum(d^2)/sum(d)-1)
    BH <- (r^2-1)*I-r*A+D
    #rho <- partial_eigen(BH,15)$values#eigs_sym(BH,15,which="LM",sigma=0)$values
    rho <- sort(eigs_sym(BH,15,which="SA")$values)
    diff <- rho[2:K.max]-5*rho[1:(K.max-1)]
    #if(rho[1]>5*rho[2]) return(TRUE)
    return(list(K=max(which(diff>0)),values=rho))
}









NCV.select <- function(A,max.K,cv=3){
dc.avg.se <- dc.avg.log <- avg.se <- avg.log <- rep(0,max.K)
dc.avg.se[1] <- dc.avg.log[1] <- avg.se[1] <- avg.log[1] <- Inf
dc.dev.mat <- dc.l2.mat <- sbm.dev.mat <- sbm.l2.mat <- matrix(0,cv,max.K)
n <- nrow(A)
sample.index <- sample.int(n)
max.fold.num <- ceiling(n/cv)
fold.index <- rep(1:cv,each=max.fold.num)[1:n]

cv.index <- fold.index[sample.index]
for(KK in 1:max.K){
    dc.l2 <- l2 <- dc.log.like <- log.like <- rep(0,cv)
    print(paste("Start",KK))
    for(k in 1:cv){
        holdout.index <- which(cv.index==k)
        train.index <- which(cv.index!=k)
            tmp.eval <- cv.evaluate(A,train.index,holdout.index,KK)
            tmp.eval.dc <- cv.evaluate.DC(A,train.index,holdout.index,KK)

        log.like[k] <- tmp.eval$loglike

        sbm.l2.mat[k,KK] <- l2[k] <- tmp.eval$l2
         sbm.dev.mat[k,KK] <- log.like[k] <- tmp.eval$loglike

        dc.l2.mat[k,KK] <- dc.l2[k] <- tmp.eval.dc$l2
        dc.dev.mat[k,KK] <- dc.log.like[k] <- tmp.eval.dc$loglike
    }

    avg.se[KK] <- mean(l2)
    avg.log[KK] <- mean(log.like)
    dc.avg.se[KK] <- mean(dc.l2)
    dc.avg.log[KK] <- mean(dc.log.like)


    print(paste("Finish ",KK,"....",sep=""))
}
   if(min(avg.log)>min(dc.avg.log)){
       dev.model <- paste("DCSBM",which.min(dc.avg.log),sep="-")
   }else{
       dev.model <- paste("SBM",which.min(avg.log),sep="-")
   }
   if(min(avg.se)>min(dc.avg.se)){
       l2.model <- paste("DCSBM",which.min(dc.avg.se),sep="-")
   }else{
       l2.model <- paste("SBM",which.min(avg.se),sep="-")
   }
   return(list(dev=avg.log,l2=avg.se,dc.dev=dc.avg.log,dc.l2=dc.avg.se,sbm.l2.mat=sbm.l2.mat,sbm.dev.mat=sbm.dev.mat,dc.l2.mat=dc.l2.mat,dc.dev.mat=dc.dev.mat,l2.model=l2.model,dev.model=dev.model))
}




cv.evaluate <- function(A,train.index,holdout.index,K){
    n <- nrow(A)
    A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
    n.holdout <- length(holdout.index)
    n.train <- n-n.holdout
    A1 <- A.new[1:n.train,]
    A1.svd <- irlba(A1+0.001,nu=K,nv=K)
    #A1.svd <- irlba(A1,nu=K,nv=K)

    if(K==1){
      A0 <- A1[1:n.train,1:n.train]
      pb <- sum(A0)/n.train^2
      if(pb < 1e-6) pb <- 1e-6
      if(pb > 1- 1e-6) pb <- 1-1e-6
      A.2 <- A.new[(n.train+1):n,(n.train+1):n]
      sum.index <- lower.tri(A.2)
      loglike <- -sum(A.2[sum.index]*log(pb)) - sum((1-A.2[sum.index])*log(1-pb))
      l2 <- sum((A.2[sum.index]-pb)^2)
      return(list(loglike=loglike,l2=l2))
    }

    V <- A1.svd$v
    km <- kmeans(V,centers=K,nstart=30,iter.max=30)

    degrees <- colSums(A1)
    no.edge <- sum(degrees==0)


    B <- matrix(0,K,K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            N.1j <- intersect(1:n.train,which(km$cluster==j))
            N.2j <- intersect((n.train+1):n,which(km$cluster==j))
            B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(length(N.1i)*length(N.1j)+length(N.1j)*length(N.2i)+length(N.1i)*length(N.2j)+1)
            #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
        }
    }
    B <- B+t(B)
    Theta <- matrix(0,n,K)
    for(i in 1:K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(length(N.1i)*(length(N.1i)-1)/2+length(N.1i)*length(N.2i)+1)
            Theta[which(km$cluster==i),i] <- 1

    }
    P.hat.holdout <- Theta[(n.train+1):n,]%*%B%*%t(Theta[(n.train+1):n,])
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))

    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2))



}






cv.evaluate.DC <- function(A,train.index,holdout.index,K){
    n <- nrow(A)
    A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
    n.holdout <- length(holdout.index)
    n.train <- n-n.holdout
    A1 <- A.new[1:n.train,]
    A1.svd <- irlba(A1,nu=K,nv=K)
    V <- A1.svd$v[,1:K]
    if(K==1) {V.norms <- abs(V)}else{
    V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    }
    iso.index <- which(V.norms==0)
    Psi <- V.norms
    Psi <- Psi / max(V.norms)
    Psi.outer <- outer(Psi,Psi)
    inv.V.norms <- 1/V.norms
    inv.V.norms[iso.index] <- 1
    V.normalized <- diag(inv.V.norms)%*%V

    if(K==1){
            N.1i <- 1:n.train
            N.2i <- (n.train+1):n
            pb <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer))+1)


    P.hat.holdout <-  diag(Psi[(n.train+1):n])%*%matrix(1,ncol=(n-n.train),nrow=(n-n.train))%*%diag(Psi[(n.train+1):n])*pb
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))

    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2))
    }


    km <- kmeans(V.normalized,centers=K,nstart=30,iter.max=30)

    degrees <- colSums(A1)
    no.edge <- sum(degrees==0)

    B <- matrix(0,K,K)

    for(i in 1:(K-1)){
        for(j in (i+1):K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            N.1j <- intersect(1:n.train,which(km$cluster==j))
            N.2j <- intersect((n.train+1):n,which(km$cluster==j))
            B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(sum(Psi.outer[N.1i,N.1j]) + sum(Psi.outer[N.1i,N.2j]) + sum(Psi.outer[N.1j,N.2i])+1)
            #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
        }
    }
    B <- B+t(B)
    Theta <- matrix(0,n,K)
    for(i in 1:K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer))+1)
            Theta[which(km$cluster==i),i] <- 1

    }
    tmp.imt.mat <- Theta[(n.train+1):n,]*Psi[(n.train+1):n]
    P.hat.holdout <-  tmp.imt.mat%*%B%*%t(tmp.imt.mat)
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))

    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2))

}








###### low rank matrix completion based on two methods used in Li et al, 2016



iter.SVD.core <- function(A,K,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,fast=FALSE,p.sample=1,kappa=NULL){
    if(sparse) A <- Matrix(A,sparse=TRUE)
     avg.p <- mean(as.numeric(A),na.rm=TRUE)
   if(is.null(kappa))kappa <- 2/avg.p
    cap <- 1#kappa*avg.p
    if(cap>1-tau) cap <- 1-tau
    if(fast){
        if(verbose) print("Matrix completion with fast approximation!")
        A[which(is.na(A))] <- 0
        A <- A/p.sample
        #svd.new <- svd(A,nu=K,nv=K)
        svd.new <- irlba(A,nu=K,nv=K)
        if(K==1){ A.new <- svd.new$d[1]*matrix(svd.new$u,ncol=1)%*%t(matrix(svd.new$v,ncol=1))}else{
        A.new <- svd.new$u%*%(t(svd.new$v)*svd.new$d[1:K])}
        A.new[A.new < 0+tau] <- 0+tau
        A.new[A.new >cap] <- cap
        return(list(iter=NA,SVD=svd.new,A=A.new,err.seq=NA))
    }
    #if(sparse) A <- Matrix(A,sparse=TRUE)
    Omega <- which(is.na(A))
    A.col.means <- colMeans(A,na.rm=TRUE)
    #avg.p <- mean(as.numeric(A),na.rm=TRUE)
    A.impute <- A
    n <- nrow(A)
    p <- ncol(A)
    if(is.null(init)){
        A.impute[Omega] <- runif(n=length(Omega))
        init.SVD <- irlba(A.impute,nu=K,nv=K)
    }else{
        init.SVD <- init
    }
    #print(init.SVD$u)
    ## init.SVD <- irlba(A,nu=K,nv=K) ## if you are working on large problems
    if(K==1){U.old <- V.old <- matrix(0,n,1)}else{
    if(K==2){
    U.old <- matrix(init.SVD$u[,1:(K-1)],ncol=K-1)*(sqrt(init.SVD$d[1:(K-1)]))
    V.old <- matrix(init.SVD$v[,1:(K-1)],ncol=K-1)*(sqrt(init.SVD$d[1:(K-1)]))
    }else{
    #print(init.SVD$u)
    U.old <- matrix(init.SVD$u[,1:(K-1)],ncol=K-1)%*%diag(sqrt(init.SVD$d[1:(K-1)]))
    V.old <- matrix(init.SVD$v[,1:(K-1)],ncol=K-1)%*%diag(sqrt(init.SVD$d[1:(K-1)]))
    }
    }
    A.old <- U.old %*% t(V.old)
    R <- A - A.old

    R[Omega] <- 0
    if(verbose) print(norm(R))
    A.impute <- A
    A.impute[Omega] <- A.old[Omega]
    ### begin iteration
    flag <- 0
    iter <- 0
    err.seq <- norm(R,"F")
    shrink <- 0
    while((iter < max.iter) && (flag != 1)){
        #print(iter)
        iter <- iter + 1
        svd.new <- irlba(A.impute,nu=K,nv=K)
        if(K==1){ A.new <- svd.new$d[1]*matrix(svd.new$u,ncol=1)%*%t(matrix(svd.new$v,ncol=1))}else{
        A.new <- svd.new$u%*%diag(svd.new$d)%*%t(svd.new$v)}
        A.new[A.new < 0+tau] <- 0+tau
        A.new[A.new >cap] <- cap
        A.impute[Omega] <- A.new[Omega]
        A.old <- A.new
        R <- A.impute - A.new
        err <- norm(R,"F")
        if(verbose) print(err)
        err.seq <- c(err.seq,err)
        if(abs(err.seq[iter+1]-err.seq[iter])<tol*err.seq[iter]) flag <- 1
     }
    #print(iter)
    return(list(iter=iter,SVD=svd.new,A=A.new,err.seq=err.seq))
}




## more time efficient (but less memory efficient) version of previous algorithm, if one uses the fast algorithm
iter.SVD.core.fast.all <- function(A,Kmax,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,p.sample=1){
    if(sparse) A <- Matrix(A,sparse=TRUE)
     avg.p <- mean(as.numeric(A),na.rm=TRUE)
    cap <- 1#kappa*avg.p
        A[which(is.na(A))] <- 0
        A <- A/p.sample
        #svd.new <- svd(A,nu=K,nv=K)
     #print("begin SVD")
        svd.new <- irlba(A,nu=Kmax,nv=Kmax)
    #print("end SVD")
        result <- list()
        for(K in 1:Kmax){
            print(K)
        if(K==1){
            A.new <- svd.new$d[1]*matrix(svd.new$u[,1],ncol=1)%*%t(matrix(svd.new$v[,1],ncol=1))
              }else{
        A.new <- A.new + svd.new$d[K]*matrix(svd.new$u[,K],ncol=1)%*%t(matrix(svd.new$v[,K],ncol=1))
        }
        A.new.thr <- A.new
        A.new.thr[A.new < 0+tau] <- 0+tau
        A.new.thr[A.new >cap] <- cap

        tmp.SVD <- list(u=svd.new$u[,1:K],v=svd.new$v[,1:K],d=svd.new$d[1:K])
        result[[K]] <- list(iter=NA,SVD=tmp.SVD,A=A.new,err.seq=NA,A.thr=A.new.thr)
        }
        return(result)

}





### ECV method for block model selection, from Li et al. 2016







ECV.block <- function(A,max.K,cv=NULL,B=3,holdout.p=0.1,tau=0,dc.est=2,kappa=NULL){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    holdout.index.list <- list()
    if(is.null(cv)){
        holdout.n <- floor(holdout.p*edge.n)

        for(j in 1:B){
            holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
        }
    }else{
        sample.index <- sample.int(edge.n)
        max.fold.num <- ceiling(edge.n/cv)
        fold.index <- rep(1:cv,each=max.fold.num)[edge.n]
        cv.index <- fold.index[sample.index]
        B <- cv
        for(j in 1:B){
            holdout.index.list[[j]] <- which(cv.index==j)
        }
    }
    #print(fast)
    result <- lapply(holdout.index.list,holdout.evaluation.fast.all,A=A,max.K=max.K,tau=tau,dc.est=dc.est,p.sample=1-holdout.p,kappa=kappa)
    dc.block.err.mat <- dc.loglike.mat <- bin.dev.mat <- roc.auc.mat <- impute.err.mat <- block.err.mat <- loglike.mat <- matrix(0,nrow=B,ncol=max.K)
    no.edge.seq <- rep(0,B)
    Omega.list <- A.list <- Imputed.A.list <- list()
    for(b in 1:B){
        impute.err.mat[b,] <- result[[b]]$impute.sq.err
        block.err.mat[b,] <- result[[b]]$block.sq.err
        loglike.mat[b,] <- result[[b]]$loglike
        roc.auc.mat[b,] <- result[[b]]$roc.auc
        bin.dev.mat[b,] <- result[[b]]$bin.dev
        no.edge.seq[b] <- result[[b]]$no.edge
        dc.block.err.mat[b,] <- result[[b]]$dc.block.sq.err
        dc.loglike.mat[b,] <- result[[b]]$dc.loglike

    }


    output <- list(impute.err=colMeans(impute.err.mat),l2=colMeans(block.err.mat),dev=colSums(loglike.mat),auc=colMeans(roc.auc.mat),dc.l2=colMeans(dc.block.err.mat),dc.dev=colSums(dc.loglike.mat),sse=colMeans(impute.err.mat),auc.mat=roc.auc.mat,dev.mat=loglike.mat,l2.mat=block.err.mat,SSE.mat=impute.err.mat,dc.dev.mat=dc.loglike.mat,dc.l2.mat=dc.block.err.mat)

   if(min(output$dev)>min(output$dc.dev)){
       dev.model <- paste("DCSBM",which.min(output$dc.dev),sep="-")
   }else{
       dev.model <- paste("SBM",which.min(output$dev),sep="-")
   }
   if(min(output$l2)>min(output$dc.l2)){
       l2.model <- paste("DCSBM",which.min(output$dc.l2),sep="-")
   }else{
       l2.model <- paste("SBM",which.min(output$l2),sep="-")
   }
  output$l2.model <- l2.model
  output$dev.model <- dev.model

    return(output)
}




holdout.evaluation.fast.all <- function(holdout.index,A,max.K,tau=0,dc.est=1,p.sample=1,kappa=NULL){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    non.miss <- which(!is.na(A.new))
    #A.new[non.miss] <- A.new[non.miss] + 0.5
    SVD.result <- iter.SVD.core.fast.all(A.new,max.K,p.sample=p.sample)
    dc.block.sq.err <-  dc.loglike <- roc.auc <- bin.dev <- block.sq.err <- impute.sq.err <- loglike <- rep(0,max.K)

    for(k in 1:max.K){
        print(k)
        #print(fast)
        tmp.est <- SVD.result[[k]]
        A.approx <- tmp.est$A.thr
        impute.sq.err[k] <- sum((A.approx[Omega]-A[Omega])^2)
        response <- A[edge.index[holdout.index]]#A[Omega]
        predictors <- A.approx[edge.index[holdout.index]]#A.approx[Omega]
        print("AUC claculation")
        #print(system.time(tmp.roc <- pROC::roc(response=response,predictor=predictors)))
        #print(length(unique(predictors)))
        aa <- roc(predictions=predictors,labels=factor(response))
        #tmp.roc.smooth <- smooth(tmp.roc,method="binormal")
        roc.auc[k] <- auc(aa)#as.numeric(tmp.roc$auc)
        #print(tmp.roc$auc)
        #print(auc(aa))
        #roc.auc[k] <- as.numeric(tmp.roc.smooth$auc)
        trunc.predictors <- predictors
        trunc.predictors[predictors>(1-1e-6)] <- 1-1e-6
        trunc.predictors[predictors<1e-6] <- 1e-6
        bin.dev[k] <- sum((response-trunc.predictors)^2)#-sum(response*log(trunc.predictors)) - sum((1-response)*log(1-trunc.predictors))
        if(k==1){
            pb <- (sum(A.new,na.rm=TRUE)+1)/(sum(!is.na(A.new)) -sum(!is.na(diag(A.new)))+1)
            if(pb < 1e-6) pb <- 1e-6
            if(pb > 1-1e-6) pb <- 1-1e-6
            A.Omega <- A[Omega]
            block.sq.err[k] <- sum((pb-A[Omega])^2)
            loglike[k] <- -sum(A.Omega*log(pb)) - sum((1-A.Omega)*log(1-pb))

        }

        #U.approx <- eigen(A.approx)$vectors[,1:k]
        print("SBM calculation")
        #print(k)
        #print(dim(tmp.est$SVD$v))
        ptm <- proc.time()
        if(k==1) {U.approx <- matrix(tmp.est$SVD$v,ncol=k)}else{
            U.approx <- tmp.est$SVD$v[,1:k]
            if(tau>0){
            A.approx <- A.approx + tau*mean(colSums(A.approx))/n
            d.approx <- colSums(A.approx)
            L.approx <- diag(1/sqrt(d.approx))%*%A.approx%*%diag(1/sqrt(d.approx))
            A.approx.svd <- irlba(L.approx,nu=k,nv=k)
            U.approx <- A.approx.svd$v[,1:k]
            }
        }

        km <- kmeans(U.approx,centers=k,nstart=30,iter.max=30)
        B <- matrix(0,k,k)
        Theta <- matrix(0,n,k)
        for(i in 1:k){
            for(j in i:k){
                N.i <- which(km$cluster==i)
                N.j <- which(km$cluster==j)
                if(i!=j){
                    B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j]))+1)
                } else{
                    #print(max(N.i))
                    #print(max(N.j))
                    #print(dim(A.new))
                   B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j])) -sum(!is.na(diag(A.new[N.i,N.j])))+1)
                }

            }
            Theta[N.i,i] <- 1
        }
        P.hat <- Theta%*%B%*%t(Theta)
        diag(P.hat) <- 0
        block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
        P.hat.Omega <- P.hat[Omega]
        A.Omega <- A[Omega]
        P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
        P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
        loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
        print(proc.time() - ptm)
#### Degree correct model
        V <- U.approx
        print("DCSBM calculation")
        ptm <- proc.time()
        #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        if(k==1) {V.norms <- as.numeric(abs(V))}else{
            V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        }

        iso.index <- which(V.norms==0)
        Psi <- V.norms
        Psi <- Psi / max(V.norms)
        inv.V.norms <- 1/V.norms
        inv.V.norms[iso.index] <- 1

        V.normalized <- diag(as.numeric(inv.V.norms))%*%V

        #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        #Psi <- V.norms
        #Psi <- Psi / max(V.norms)
        #V.normalized <- diag(1/V.norms)%*%V
        #Psi.outer <- outer(Psi,Psi)
        if(k==1){
        if(dc.est>1){
            B <- sum(A.new,na.rm=TRUE)+0.01

            partial.d <- colSums(A.new,na.rm=TRUE)
            partial.gd <- B
            phi <- rep(0,n)
            B.g <- partial.gd
            phi <- as.numeric(partial.d/B.g)
            B <- B/p.sample
            P.hat <- t(t(matrix(B,n,n)*phi)*phi)
            #P.hat <- diag(phi)%*%matrix(B,n,n)%*%diag(phi)
            diag(P.hat) <- 0
        }
            dc.block.sq.err[k] <- sum((pb-A[Omega])^2)
            P.hat.Omega <- P.hat[Omega]
            A.Omega <- A[Omega]
            P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
            P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6

            dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))


        }else{
        km <- kmeans(V.normalized,centers=k,nstart=30,iter.max=30)
        if(dc.est>1){
            B <- matrix(0,k,k)
            Theta <- matrix(0,n,k)
            for(i in 1:k){
                for(j in 1:k){
                N.i <- which(km$cluster==i)
                N.j <- which(km$cluster==j)
                B[i,j] <- sum(A.new[N.i,N.j],na.rm=TRUE)+0.01
                }
                Theta[N.i,i] <- 1
            }
            Theta <- Matrix(Theta,sparse=TRUE)
            partial.d <- colSums(A.new,na.rm=TRUE)
            partial.gd <- colSums(B)
            phi <- rep(0,n)
            B.g <- Theta%*%partial.gd
            phi <- as.numeric(partial.d/B.g)
            B <- B/p.sample
            tmp.int.mat <- Theta*phi
            P.hat <-as.matrix(tmp.int.mat%*%B%*%t(tmp.int.mat))
            #P.hat <- diag(phi)%*%Theta%*%B%*%t(Theta)%*%diag(phi)
            diag(P.hat) <- 0
        }
        dc.block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
        P.hat.Omega <- P.hat[Omega]
        A.Omega <- A[Omega]
        P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
        P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
        dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
    }
        print(proc.time() - ptm)



    }
    return(list(impute.sq.err=impute.sq.err,block.sq.err=block.sq.err,loglike=loglike,roc.auc=roc.auc,no.edge=no.edge,dc.block.sq.err=dc.block.sq.err,dc.loglike=dc.loglike,bin.dev=bin.dev))
}










DirectedRDPG <- function(n,K,avg.d=NULL){
    Z1 <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    Z2 <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    S <- Z1%*%t(Z2)
    P <- S/max(S)
    if(!is.null(avg.d)){
        P <- P*avg.d/(mean(rowSums(P)))
    }

    A <- matrix(0,n,n)
    R <- matrix(runif(n^2),n,n)
    A[R<P] <- 1
    return(list(A=A,P=P))
}







UndirectedRDPG <- function(n,K,avg.d=NULL){
    Z1 <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    S <- Z1%*%t(Z1)
    P <- S/max(S)
    if(!is.null(avg.d)){
        P <- P*avg.d/(mean(rowSums(P)))
    }

    upper.index <- which(upper.tri(P))
    upper.p <- P[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1


    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(list(A=A,P=P))
}

RDPG.Gen <- function(n,K,directed=TRUE,avg.d=NULL){
    if(directed){
        return(DirectedRDPG(n,K,avg.d))
    }else{
        return(UndirectedRDPG(n,K,avg.d))
    }
}








#### use imputation error to estimate the rank of an undirected network. Note that if select weighted = TRUE, AUC will not be calculated and the computation will be much faster.

ECV.undirected.Rank <- function(A,max.K,B=3,holdout.p=0.1,weighted=TRUE){
    n <- nrow(A)
    #edge.index <- 1:n^2
    #edge.n <- length(edge.index)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)

    holdout.index.list <- list()

    holdout.n <- floor(holdout.p*edge.n)

    for(j in 1:B){
        holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
    }
    result <- lapply(holdout.index.list,missing.undirected.Rank.fast.all,A=A,max.K=max.K,p.sample=1-holdout.p,weighted)
    sse.mat <- roc.auc.mat <- matrix(0,nrow=B,ncol=max.K)

    for(b in 1:B){
        roc.auc.mat[b,] <- result[[b]]$roc.auc
        sse.mat[b,] <- result[[b]]$sse
    }
    if(!weighted){
    auc.seq <- colMeans(roc.auc.mat)
    }else{
        auc.seq <- rep(NA,max.K)
    }
    sse.seq <- colMeans(sse.mat)
    return(list(sse.rank=which.min(sse.seq),auc.rank=which.max(auc.seq),auc=auc.seq,sse=sse.seq))
}





missing.undirected.Rank.fast.all <- function(holdout.index,A,max.K,p.sample=1,weighted=TRUE){
    n <- nrow(A)
    #A.new <- A
    #A.new[holdout.index] <- NA
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    diag(A.new) <- diag(A)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    sse <- roc.auc <- rep(0,max.K)
    SVD.result <- iter.SVD.core.fast.all(A.new,max.K,p.sample=p.sample)
    for(k in 1:max.K){
        print(k)
        tmp.est <- SVD.result[[k]]
        #if(k==1){
        #A.approx <- matrix(tmp.est$SVD$u,ncol=1)%*%t(matrix(tmp.est$SVD$v,ncol=1))*tmp.est$SVD$d[1]
        #}else{
         #   A.approx <- tmp.est$SVD$u%*%t(tmp.est$SVD$v*tmp.est$SVD$d)
        #}
        A.approx <- tmp.est$A
        response <- A[Omega]
        predictors <- A.approx[Omega]
        if(!weighted){
        aa <- roc(predictions=predictors,labels=factor(response))
        roc.auc[k] <- auc(aa)
        }
        sse[k] <- mean((response-predictors)^2)
        imputed.A[[k]] <- A.approx
    }
    return(list(imputed.A=imputed.A,Omega=Omega,sse=sse,roc.auc=roc.auc))
}


ECV.directed.Rank <- function(A,max.K,B=3,holdout.p=0.1,weighted=TRUE){
    n <- nrow(A)
    edge.index <- 1:n^2
    edge.n <- length(edge.index)
    #edge.index <- which(upper.tri(A))
    #edge.n <- length(edge.index)

    holdout.index.list <- list()

    holdout.n <- floor(holdout.p*edge.n)

    for(j in 1:B){
        holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
    }
    result <- lapply(holdout.index.list,missing.directed.Rank.fast.all,A=A,max.K=max.K,p.sample=1-holdout.p,weighted=weighted)
    sse.mat <- roc.auc.mat <- matrix(0,nrow=B,ncol=max.K)

    for(b in 1:B){
        roc.auc.mat[b,] <- result[[b]]$roc.auc
        sse.mat[b,] <- result[[b]]$sse
    }

    if(!weighted){
    auc.seq <- colMeans(roc.auc.mat)
    }else{
        auc.seq <- rep(NA,max.K)
    }
    sse.seq <- colMeans(sse.mat)
    return(list(sse.rank=which.min(sse.seq),auc.rank=which.max(auc.seq),auc=auc.seq,sse=sse.seq))
    #return(list(rank=which.max(auc.seq),auc=auc.seq))
}





missing.directed.Rank.fast.all <- function(holdout.index,A,max.K,p.sample=NULL,weighted=TRUE){
    n <- nrow(A)
    A.new <- A
    A.new[holdout.index] <- NA
    if(is.null(p.sample)){
        p.sample <- 1-length(holdout.index)/n^2
    }
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    sse <- roc.auc <- rep(0,max.K)
    SVD.result <- iter.SVD.core.fast.all(A.new,max.K,p.sample=p.sample)

    for(k in 1:max.K){
        print(k)
        tmp.est <- SVD.result[[k]]
        A.approx <- tmp.est$A.thr
        response <- A[Omega]
        predictors <- A.approx[Omega]
        if(!weighted){
        aa <- roc(predictions=predictors,labels=factor(response))
        roc.auc[k] <- auc(aa)
        }
        sse[k] <- mean((response-predictors)^2)
        imputed.A[[k]] <- A.approx
    }
    return(list(roc.auc=roc.auc,imputed.A=imputed.A,Omega=Omega,sse=sse))
}

ECV.Rank <- function(A,max.K,B=3,holdout.p=0.1,weighted=TRUE,mode="directed"){
    if(mode=="directed"){
        return(ECV.directed.Rank(A=A,max.K=max.K,B=B,holdout.p=holdout.p,weighted=weighted))
    }else{
        return(ECV.undirected.Rank(A=A,max.K=max.K,B=B,holdout.p=holdout.p,weighted=weighted))
    }
}



##### neighborhood smoothing estimation of graphon model from Zhang et al. 2015


nSmooth <- function(A,h=NULL){
    n <- nrow(A)
    A2 <- A%*%A/n
    if(is.null(h)){
        h <- sqrt(log(n)/n)
    }
    Kmat <- D <- matrix(0,n,n)
    for(k in 1:n){
        tmp <- abs(A2 - A2[,k])
        tmp.d <- apply(tmp,2,max)
        Kmat[k,which(tmp.d < quantile(tmp.d,h))] <- 1
    }
    Kmat <- Kmat/(rowSums(Kmat)+1e-10)
    Phat <- Kmat%*%A
    Phat <- (Phat+t(Phat))/2
    return(Phat)

}


ECV.nSmooth.lowrank <- function(A,h.seq,K,cv=NULL,B=3,holdout.p=0.1){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    holdout.index.list <- list()
    if(is.null(cv)){
        holdout.n <- floor(holdout.p*edge.n)

        for(j in 1:B){
            holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
        }
    }else{
        sample.index <- sample.int(edge.n)
        max.fold.num <- ceiling(edge.n/cv)
        fold.index <- rep(1:cv,each=max.fold.num)[edge.n]
        cv.index <- fold.index[sample.index]
        B <- cv
        for(j in 1:B){
            holdout.index.list[[j]] <- which(cv.index==j)
        }
    }
    #print(fast)
    result <- lapply(holdout.index.list,holdout.evaluation.graphon.lowrank,A=A,h.seq=h.seq,K=K,p.sample=1-holdout.p)
    est.err.mat <-  matrix(0,nrow=B,ncol=length(h.seq))
    for(b in 1:B){
        est.err.mat[b,] <- result[[b]]
    }

    return(list(err=colMeans(est.err.mat),min.index=which.min(colMeans(est.err.mat))))
}



holdout.evaluation.graphon.lowrank <- function(holdout.index,A,h.seq,p.sample=1,K){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    non.miss <- which(!is.na(A.new))

    #A.new[Omega] <- 0

    SVD <- iter.SVD.core(A=A.new,K=K,p.sample=p.sample)
    A.approx <- SVD$A
   err.seq <- rep(0,length(h.seq))

   for(k in 1:length(h.seq)){
       What <- nSmooth(A.approx,h=h.seq[k])
       err.seq[k] <- mean((A[Omega]-What[Omega])^2)
   }



    return(err.seq)
}

