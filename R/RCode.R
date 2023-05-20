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
    P <- P0

    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    M[cbind(1:n,membership)] <- 1
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
    A.bar <- t(t(A.bar*node.degree)*node.degree)
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



reg.SP <- function(A,K,tau=1,lap=FALSE,nstart=30,iter.max=100){
    avg.d <- mean(colSums(A))
    A.tau <- A + tau*avg.d/nrow(A)
    if(!lap){SVD <- irlba(A.tau,nu=K,nv=K)}else{
         d.tau <- colSums(A.tau)
         L.tau <- diag(1/sqrt(d.tau))%*%A.tau%*%diag(1/sqrt(d.tau))
         #SVD <- svd(L.tau,nu=K,nv=K)
         SVD <- irlba(L.tau,nu=K,nv=K)
    }
    km <- kmeans(SVD$v[,1:K],centers=K,nstart=nstart,iter.max=iter.max)#,algorithm="Lloyd")
    return(list(cluster=km$cluster,loss=km$tot.withinss))
}


reg.SSP <- function(A,K,tau=1,lap=FALSE,nstart=30,iter.max=100){
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
    km <- kmeans(V.normalized,centers=K,nstart=nstart,iter.max=iter.max)#,algorithm="Lloyd")
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
                if(n.i>1){
                  B[i,i] <- sum(A[which(g==i),which(g==i)])/(n.i^2 - n.i)
                }else{
                  B[i,i] <- 0
                }
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
            model <- "DCSBM"
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
    rho <- sort(eigs_sym(BH,K.max,which="SA")$values)
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
            #print(K)
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
        #print(k)
        #print(fast)
        tmp.est <- SVD.result[[k]]
        A.approx <- tmp.est$A.thr
        impute.sq.err[k] <- sum((A.approx[Omega]-A[Omega])^2)
        response <- A[edge.index[holdout.index]]#A[Omega]
        predictors <- A.approx[edge.index[holdout.index]]#A.approx[Omega]
        #print("AUC claculation")
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
        #print("SBM calculation")
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
        #print(proc.time() - ptm)
#### Degree correct model
        V <- U.approx
        #print("DCSBM calculation")
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
        #print(proc.time() - ptm)



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
        #print(k)
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
        #print(k)
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




ConsensusClust <- function(A,K,tau=0.25,lap=TRUE){

  n <- nrow(A)
  sigma.list <- list()
  for(u in 1:n){
    Au <- A[-u,-u]
    sc <- reg.SP(A=Au,K=K,tau=tau,lap=lap)
    sigma.u <- sc$cluster
    est.u <- SBM.estimate(A=Au,g=sigma.u)
    a.u <- n*min(diag(est.u$B)+0.1/n)
    fake.B.u <- est.u$B
    diag(fake.B.u) <- -1
    b.u <- n*max(fake.B.u)
    t.u <- 0.5*log(a.u*(1-b.u/n))-0.5*log(b.u*(1-a.u/n))
    rho.u <- -0.5*(1/t.u)*(log(a.u/n*exp(-t.u)+1-a.u/n)-log(b.u/n*exp(t.u)+1-b.u/n))
    Z.u <- matrix(0,nrow=n-1,ncol=K)
    Z.u[cbind(1:(n-1),sigma.u)] <- 1
    u.comm <- matrix(A[u,-u],nrow=1)%*%Z.u
    u.size <- colSums(Z.u)
    score <- u.comm - rho.u*u.size
    tmp.sigma <- rep(0,n)
    tmp.sigma[-u] <- sigma.u
    tmp.sigma[u] <- which.max(score)
    sigma.list[[u]] <- tmp.sigma
  }
  ## consensus step
  sigma <- rep(0,n)
  sigma[1] <- sigma.list[[1]][1]
  base.sigma <- sigma.list[[1]]
  for(u in 2:n){
    K.score <- rep(0,K)
    for(k in 1:K){
      K.score[k] <- length(intersect(which(base.sigma==k),which(sigma.list[[u]]==sigma.list[[u]][u])))
    }
    sigma[u] <- which.max(K.score)
  }
  return(sigma)
}




RightSC <- function(A,K,normal=FALSE){
  SVD <- irlba(A,nv=K,nu=K)
  V <- SVD$v
  if(normal){
    Vnorm <- apply(V,1,function(x)sqrt(sum(x^2)))
    V <- V/(0.001+Vnorm)
  }
  km <- kmeans(V,centers=K,iter.max=200,nstart=50)
  return(km)
}


NSBM.estimate <- function(A,K,g,reg.bound=-Inf){

  n <- nrow(A)
  B <- matrix(0,K,K)
  diag(B) <- 1
  theta <- rep(1,n)
  Z <- matrix(0,n,K)
  Z[cbind(1:n,g)] <- 1
  N2C <- A%*%Z+1
  nK <- as.numeric(table(g))
  N2C <- t(t(N2C)/nK)
  theta <- N2C[cbind(1:n,g)]
  Y <- log(N2C)
  membership.list <- list()
  for(k in 1:K){
    membership.list[[k]] <- which(g==k)
  }
  ##### estimate B
  for(k in 1:K){
    for(l in 1:K){
      if(l!=k){
        tmp <- Y[,k]-Y[,l]
        B[k,l] <- exp(-mean(tmp[membership.list[[k]]]))
      }
    }
  }

  lambda <- rep(1,n)
  for(k in 1:K){
    tmp.mat <- -(Y - Y[,k])
    tmp.sum <- colMeans(tmp.mat[membership.list[[k]],])
    tmp <- rowSums(tmp.mat)
    tmp2 <- sum(tmp.sum)
    lambda[membership.list[[k]]] <- tmp[membership.list[[k]]]/tmp2
  }

  #reg.bound <- -2 ### some regularity lower bound for lambda,for numerical stability
  if(min(lambda)<reg.bound){
    for(k in 1:K){
      index <- which(g==k)
      tmp.lambda <- lambda[index]
      incons.index <- which(tmp.lambda< reg.bound)
      if(length(incons.index)>0){
        cons.index <- which(tmp.lambda>= reg.bound)
        tmp.lambda[incons.index] <- -reg.bound
        tmp.lambda[cons.index] <- tmp.lambda[cons.index]*((length(index)-reg.bound*length(incons.index))/sum(tmp.lambda[cons.index]))
        lambda[index] <- tmp.lambda
      }
    }
  }


  P.tilde <- Z%*%B%*%t(Z)
  for(i in 1:n){
    P.tilde[i,] <- (P.tilde[i,]^lambda[i])*theta[i]
  }
  return(list(B=B,lambda=lambda,theta=theta,P.tilde=P.tilde,g=g))

}







NSBM.Gen <- function(n,K,avg.d,beta,theta.low=0.1,theta.p=0.2,lambda.scale=0.2,lambda.exp=FALSE){
  membership <- sample(K,size=n,replace=TRUE)
  if(length(beta)==1){
    B <- matrix(beta,K,K)
    diag(B) <- 1
  }else{
    B <- beta
  }
  Z <- matrix(0,n,K)
  Z[cbind(1:n,membership)] <- 1
  P <- Z%*%B%*%t(Z)
  if(!lambda.exp){
    lambda <- rnorm(n)*lambda.scale + 1
    while(sum(lambda<0)>0.05){
      lambda <- rnorm(n)*lambda.scale + 1
    }
  }else{
    log.lambda <- lambda.scale*2*runif(n)-lambda.scale
    lambda <- exp(log.lambda)
  }
  theta <- rep(1,n)
  theta[sample(n,size=ceiling(n*theta.p))] <- theta.low
  for(k in 1:K){
    gk <- which(membership==k)
    lambda[gk] <- lambda[gk]-mean(lambda[gk])+1
  }
  P.tilde <- P
  for(i in 1:n){
    P.tilde[i,] <- ((P[i,]^lambda[i])*theta[i])
  }
  current.avg.d <- mean(rowSums(P.tilde))
  P.tilde <- P.tilde/current.avg.d*avg.d
  theta <- theta/current.avg.d*avg.d
  A <- matrix(0,n,n)
  A[matrix(runif(n^2),n,n) < P.tilde] <- 1
  return(list(A=A,P=P,P.tilde=P.tilde,B=B,theta=theta,lambda=lambda,g=membership))
}




LSM.PGD <- function(A,k,step.size=0.3,niter=500,trace=0){
  N <- nrow(A)
  ones = rep(1,N)
  M = matrix(1, N, N)
  Jmat <- diag(rep(1,N)) - M/N

  P.tilde <- USVT(A)
  P.tilde[P.tilde>(1-1e-5)] <- (1-1e-5)
  P.tilde[P.tilde< 1e-5] <- 1e-5

  Theta.tilde <- logit(P.tilde)

  alpha_0 <- solve(N*diag(rep(1,N))+M,rowSums(Theta.tilde))

  G <- Jmat%*%(Theta.tilde - outer(alpha_0,alpha_0,"+"))%*%Jmat

  eig <- eigs_sym(A=G,k = k)
  eig$values[eig$values<=0] <- 0

  Z_0 <- t(t(eig$vectors[,1:k])*sqrt(eig$values[1:k]))
  obj <- NULL
  step.size.z <- step.size/norm(Z_0,"2")^2
  step.size.alpha <- step.size/(2*N)
  for(i in 1:niter){
    Theta.hat <- alpha_0 %*% t(rep(1,N)) + rep(1, N) %*% t(alpha_0) + Z_0 %*% t(Z_0)
    Phat <- sigmoid(Theta.hat)
    tmp.obj <- (sum(A*log(Phat)) + sum((1-A)*log(1-Phat)) - sum(diag(log(1-Phat))))/2
    if(trace>0){
      print(tmp.obj)
    }
    obj <- c(obj,tmp.obj)
    Z <- Z_0 + 2*step.size.z*(A-Phat)%*%Z_0
    alpha <- alpha_0 + 2*step.size.alpha*(A-Phat)%*%matrix(rep(1,N))
    Z <- Jmat%*%Z

    Z_0 <- Z
    alpha_0 <- alpha
  }

  Theta.hat <- alpha_0 %*% t(rep(1,N)) + rep(1, N) %*% t(alpha_0) + Z_0 %*% t(Z_0)
  Phat <- sigmoid(Theta.hat)
  tmp.obj <- (sum(A*log(Phat)) + sum((1-A)*log(1-Phat)) - sum(diag(log(1-Phat))))/2
  obj <- c(obj,tmp.obj)
  return(list(Z=Z,alpha=alpha,Phat=Phat,obj=obj))

}


###  example:
# dt <- RDPG.Gen(n=600,K=2,directed=TRUE)
#
# A <- dt$A
#
# fit <- LSM.PGD(A,2)




USVT <- function(A){
  n <- nrow(A)
  K <- ceiling(n^(1/3))
  SVD <- irlba(A,nv=K,nu=K)
  Ahat <- SVD$u%*%(t(SVD$v)*SVD$d)
  Ahat[Ahat>1] <- 1
  Ahat[Ahat<0] <- 0
  return(Ahat)
}


###  example:
# dt <- RDPG.Gen(n=600,K=2,directed=TRUE)
#
# A <- dt$A
#
# fit <- USVT(A)






sbm.fit = function(A,p,max.K=2){
  ##fitting the network using SBM
  ##index - the set of entries for testing
  n <- nrow(A)

  A.partial <- A
  SVD <- irlba((A.partial+mean(colSums(A.partial))*0.05/n)/(1-p),nv=max.K,nu=max.K)
  SBM.Phat.list <- list()

  for(k in 1:max.K){
    km <- kmeans(SVD$u[,1:k],centers=k, nstart = 50, iter.max = 50)
    SBM.Phat <- SBM.estimate(A.partial,g=km$cluster)$Phat/(1-p)
    if(sum(is.na(SBM.Phat))>0){
      SBM.Phat <- matrix(0,n,n)
    }
    SBM.Phat.list[[k]] <- SBM.Phat
  }
  return(SBM.Phat.list)
}
dcsbm.fit = function(A,p,max.K=2){
  ##fitting the network using DCSBM
  n <- nrow(A)
  A.partial <- A

  SVD <- irlba((A.partial+mean(colSums(A.partial))*0.05/n)/(1-p),nv=max.K,nu=max.K)
  DCSBM.Phat.list <- list()

  for(k in 1:max.K){
    V <- matrix(SVD$v[, 1:k],ncol=k)
    V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
    V.normalized <- diag(1/V.norm) %*% V
    km <- kmeans(V.normalized, centers = k, nstart = 50, iter.max = 50)
    DCSBM.Phat <- DCSBM.estimate(A.partial,g=km$cluster)$Phat/(1-p)
    if(sum(is.na(DCSBM.Phat))>0){
      DCSBM.Phat <- matrix(0,n,n)
    }
    DCSBM.Phat.list[[k]] <- DCSBM.Phat
  }
  return(DCSBM.Phat.list)
}



network.mixing <- function(A,index=NULL,rho = 0.1,max.K=15,dcsbm=TRUE, usvt=TRUE,ns=FALSE,lsm=FALSE,lsm.k=4,trace=FALSE){
  n <- nrow(A)
  if(is.null(index)){
    upper.index <- which(upper.tri(A))
    n.upper <- length(upper.index)
    sample.index <- sample(upper.index,size=ceiling(rho*n.upper))
    tmp.A <- matrix(0,n,n)
    tmp.A[sample.index] <- NA
    tmp.A <- tmp.A+t(tmp.A)
    index <- which(is.na(tmp.A))
  }
  A.partial <- A
  A.partial[index] <- 0
  Y <- A[index]
  p <- length(index)/n^2

  time1<- system.time(SVD <- irlba((A.partial+mean(colSums(A.partial))*0.05/n)/(1-p),nv=max.K,nu=max.K))
  #print(time1)
  #SBM.Phat.list <- list()
  #DCSBM.Phat.list <- list()
  SBM.full.mat <- matrix(NA,nrow = length(A.partial),ncol=max.K)
  SBM.Xmat <- matrix(0,ncol=max.K,nrow=length(index))
  if(dcsbm){
  DCSBM.full.mat  <- matrix(NA,nrow = length(A.partial),ncol=max.K)
  DCSBM.Xmat <-  matrix(0,ncol=max.K,nrow=length(index))
  }
  if(trace) print("Fitting block models....")
  ptm <- proc.time()
  for(k in 1:max.K){

    ## SBM fit
    #print("SBM")
    ptm <- proc.time()
    km <- kmeans(SVD$u[,1:k],centers=k, nstart = 50, iter.max = 50)
    SBM.Phat <- SBM.estimate(A.partial,g=km$cluster)$Phat/(1-p)
    if(sum(is.na(SBM.Phat))>0){
      SBM.Phat <- matrix(0,n,n)
    }
    #SBM.Phat.list[[k]] <- SBM.Phat
    SBM.full.mat[,k] <- as.numeric(SBM.Phat)
    SBM.Xmat[,k] <- SBM.Phat[index]
    #print(proc.time() - ptm)


    if(dcsbm){
      ## DCSBM fit
      #print("DCSBM")
      #ptm <- proc.time()
      V <- matrix(SVD$v[, 1:k],ncol=k)
      V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
      V.normalized <- diag(1/V.norm) %*% V
      km <- kmeans(V.normalized, centers = k, nstart = 50, iter.max = 50)
      DCSBM.Phat <- DCSBM.estimate(A.partial,g=km$cluster)$Phat/(1-p)
      if(sum(is.na(DCSBM.Phat))>0){
        DCSBM.Phat <- matrix(0,n,n)
      }
      #DCSBM.Phat.list[[k]] <- DCSBM.Phat
      DCSBM.full.mat[,k] <- as.numeric(DCSBM.Phat)
      DCSBM.Xmat[,k] <- DCSBM.Phat[index]
    }
    #print(proc.time() - ptm)

  }
  if(dcsbm){
    model.names <- c(paste("SBM",1:max.K,sep=""),paste("DCSBM",1:max.K,sep=""))
    Xmat <- cbind(SBM.Xmat,DCSBM.Xmat)
    full.mat <- cbind(SBM.full.mat,DCSBM.full.mat)
  }else{
    model.names <- paste("SBM",1:max.K,sep="")
    Xmat <- SBM.Xmat
    full.mat <- SBM.full.mat

  }
  #### USVT fit
  if(usvt){
    if(trace) print("USVT....")
    time2 <- system.time(usvt.est <- USVT(A.partial)/(1-p))
    #print(time2)
    Xmat <- cbind(Xmat,usvt.est[index])
    full.mat <- cbind(full.mat,as.numeric(usvt.est))
    model.names <- c(model.names,"USVT")
  }

  #### Neighborhood smoothing fit
  if(ns){
    if(trace) print("Neighborhood Smoothing....")
    ns.est <- nSmooth(A.partial)/(1-p)
    Xmat <- cbind(Xmat,ns.est[index])
    full.mat <- cbind(full.mat,as.numeric(ns.est))
    model.names <- c(model.names,"NS")
  }


  #### Latent space model fit
  if(lsm){
    if(trace)  print("Latent space model ....")
    lsm.est <- LSM.PGD(as.matrix(A.partial),lsm.k,step.size=0.3,niter=50,trace=0)$Phat/(1-p)
    Xmat <- cbind(Xmat,lsm.est[index])
    full.mat <- cbind(full.mat,as.numeric(lsm.est))
    model.names <- c(model.names,"LSM")
  }


  ## exponential aggregation
  X.l2 <- colSums((Xmat-Y)^2)
  center.l2 <- X.l2 - quantile(X.l2,0.2)
  exp.weight <- exp(-center.l2)
  exp.weight <- exp.weight/sum(exp.weight)
  #exp.Yhat <- Xmat%*%matrix(exp.weight,ncol=1)
  #exp.Yhat[exp.Yhat<0] <- 0
  #exp.Yhat[exp.Yhat>1] <- 1
  exp.Phat <- matrix(full.mat%*%matrix(exp.weight,ncol=1),ncol=n,nrow=n)
  exp.Phat[exp.Phat<0] <- 0
  exp.Phat[exp.Phat>1] <- 1

  ## model selection aggregation (ECV)
  ecv.weight <- rep(0,length(exp.weight))
  ecv.weight[which.max(exp.weight)] <- 1
  #ecv.Yhat <- Xmat%*%matrix(ecv.weight,ncol=1)
  ecv.Phat <- matrix(full.mat[, which.max(exp.weight)],nrow=n,ncol=n)
  #ecv.Yhat[ecv.Yhat<0] <- 0
  #ecv.Yhat[ecv.Yhat>1] <- 1
  ecv.Phat[ecv.Phat<0] <- 0
  ecv.Phat[ecv.Phat>1] <- 1


  #print("Fit OLS")
  #ptm <- proc.time()
  ## linear aggregation
  mixing.df <- data.frame(Xmat)
  mixing.df$Y <- Y
  lm.fit <- lm(Y~.-1,data=mixing.df)
  linear.weight <- lm.fit$coefficients
  linear.weight[which(is.na(lm.fit$coefficients))] <- 0
  #linear.Yhat <- Xmat%*%matrix(linear.weight,ncol=1)
  #linear.Yhat[linear.Yhat>1] <- 1
  #linear.Yhat[linear.Yhat<0] <- 0
  linear.Phat <- matrix(full.mat%*%matrix(linear.weight,ncol=1),ncol=n,nrow=n)
  linear.Phat[linear.Phat>1] <- 1
  linear.Phat[linear.Phat<0] <- 0
  #print(proc.time()-ptm)



  #print("Fit nnls")
  #ptm <- proc.time()
  ## nonnegative least square
  nnl <- nnls(A=Xmat,b=Y)
  nnl.weight <- nnl$x
  #nnl.Yhat <- Xmat%*%matrix(nnl.weight,ncol=1)
  #nnl.Yhat[nnl.Yhat>1] <- 1
  #nnl.Yhat[nnl.Yhat<0] <- 0
  nnl.Phat <- matrix(full.mat%*%matrix(nnl.weight,ncol=1),ncol=n,nrow=n)
  nnl.Phat[nnl.Phat<0] <- 0
  nnl.Phat[nnl.Phat>1] <- 1
  #print(proc.time()-ptm)


  return(list(linear.Phat=linear.Phat,linear.weight=linear.weight,nnl.Phat=nnl.Phat,nnl.weight=nnl.weight,exp.Phat=exp.Phat,exp.weight=exp.weight,ecv.Phat=ecv.Phat,ecv.weight=ecv.weight,model.names=model.names))
}

# dt <- RDPG.Gen(n=600,K=2,directed=TRUE)
#
# A <- dt$A
#
# fit <- network.mixing(A)
# fit$model.names
# fit$nnl.weight




network.mixing.Bfold <- function(A,B=10,rho = 0.1,max.K=15,dcsbm=TRUE,usvt=TRUE,ns=FALSE,lsm=FALSE,lsm.k=4){
  n <- nrow(A)
  upper.index <- which(upper.tri(A))
  n.upper <- length(upper.index)
  rf.Phat <- NULL
  linear.Phat=linear.Phat <- nnl.Phat <- exp.Phat <- ecv.Phat <- matrix(0,n,n)
  for(b in 1:B){
    #print(paste("Fold",b))
    sample.index <- sample(upper.index,size=ceiling(rho*n.upper))
    tmp.A <- matrix(0,n,n)
    tmp.A[sample.index] <- NA
    tmp.A <- tmp.A+t(tmp.A)
    index <- which(is.na(tmp.A))
    tmp.fit <- network.mixing(A=A,index=index,max.K=max.K,rho = rho,ns=ns,dcsbm=dcsbm, usvt=usvt,lsm=lsm,lsm.k=lsm.k)
    linear.Phat <- linear.Phat + tmp.fit$linear.Phat/B
    nnl.Phat <- nnl.Phat + tmp.fit$nnl.Phat/B
    exp.Phat <- exp.Phat + tmp.fit$exp.Phat/B
    ecv.Phat <- ecv.Phat + tmp.fit$ecv.Phat/B

  }
  return(list(linear.Phat=linear.Phat,nnl.Phat=nnl.Phat,exp.Phat=exp.Phat,ecv.Phat=ecv.Phat,model.names=tmp.fit$model.names))

}




InformativeCore <- function(A,r=3){
  n <- nrow(A)
  r = max(r, 2)
  degree = rowSums(A)
  iso.idx = (degree==0)
  er.score <- config.score <- rep(0, nrow(A))

  A = A[!iso.idx, !iso.idx]

  svd.fit = irlba(A, nv=min(r+1, nrow(A)), maxit = 200)
  A.recon = svd.fit$u[,1:r,drop=FALSE]  %*% (t( svd.fit$v[,1:r,drop=FALSE] )*svd.fit$d[1:r])
  score = apply(A.recon, 1, sd) * sqrt(nrow(A.recon)-1)
  er.score[!iso.idx] <- score

  A.recon.normal = t(t(A.recon)/rowSums(A))#A.recon %*% diag( 1/rowSums(A) )
  score = apply(A.recon.normal, 1, sd) * sqrt(nrow(A.recon.normal)-1)
  config.score[!iso.idx] <- score

  phat <- (sum(A))/(n^2-n)
  er.theory.thr <- sqrt(0.99*phat*log(n))
  er.theory.core <- which(er.score > er.theory.thr)

  config.theory.thr <- sqrt(log(n))/(n*sqrt(phat^1.01))
  config.theory.core <- which(config.score > config.theory.thr)

  er.kmeans <- kmeans(er.score,centers=2,iter.max=100,nstart=50)
  er.kmeans.core <- which(er.kmeans$cluster==which.max(er.kmeans$centers))

  config.kmeans <- kmeans(config.score,centers=2,iter.max=100,nstart=50)
  config.kmeans.core <- which(config.kmeans$cluster==which.max(config.kmeans$centers))



  return(list(er.score=er.score,config.score=config.score,er.theory.core=er.theory.core,er.kmeans.core=er.kmeans.core,config.theory.core=config.theory.core,config.kmeans.core=config.kmeans.core))
}









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




## This part of the code is modified from the contribution of an anonymous reviewer.
## load code snippet from R package 'sparseFLMM' for efficient smooth covariance estimation [Cederbaum, J., Scheipl, F., & 
## Greven, S. (2018). Fast symmetric additive covariance smoothing. Computational Statistics & Data Analysis, 120, 25-41.]
#source("symmetric_smoothing.R")

smooth.oracle <- function(Us,A){
  A <- as.matrix(A)
  N <- nrow(A)
  ## transform data to long table format
  data1 = CJ(u_1=Us, u_2=Us, sorted=FALSE)
  data1$y = c(A)
  tmpu_1 = data1$u_1
  tmpu_2 = data1$u_2
  data1$u_1[tmpu_1 > tmpu_2] = data1$u_2[tmpu_1 > tmpu_2]
  data1$u_2[tmpu_1 > tmpu_2] = tmpu_1[tmpu_1 > tmpu_2]
  data1 = data1[!is.null(data1$y)]
  
  ## spline-based graphon estimation
  estGraphon = bam(
    y ~ s(u_1, u_2, k = floor(sqrt(N)), bs = "symm", m = c(1,1), xt = list(bsmargin = 'ps')),
    family = binomial,
    data = data1,
    discrete = T
  )
  
  ## calculate predictions of edge probabilities based on graphon model fit
  data1$p_hat = predict.bam(
    estGraphon,
    newdata = data1[,c(1,2)],
    type = "response"
  )
  
  
  ## generate estimated edge probability matrix
  P_hat = matrix(NA, nrow=N, ncol=N)
  
  for(i_ in (1:N)) {
    P_hat[i_,] = data1$p_hat[(((i_ - 1) * N) + 1):(i_ * N)]
  }
  return(P_hat)
}