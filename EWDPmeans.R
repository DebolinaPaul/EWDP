wt.euc.dist.sq=function(x1,x2,w){
  p=(x1-x2)^2
  p=w*p
  return(sum(p))
}
vec.euc.dist.sq=function(x1,x2,w){
  p=(x1-x2)^2
  p=w*p
  return(p)
}
ERAkmeans=function(X,lambda_w,lambda_k,tmax,tmax1=5){
  #Initialization
  n=dim(X)[1]
  p=dim(X)[2]
  k=1
  weight=numeric(p)
  for(j in 1:p){
    weight[j]=1/p
  }
  weight=weight/sum(weight)
  label=numeric(n)
  dist=numeric(k)
  t=0
  D=numeric(p)
  Z=colMeans(X)
  # Z=X[sample(n,1),]
  dist1=numeric(k)
  #repeat the process tmax times
  repeat{
    t=t+1
    if(is.vector(Z)==TRUE){
      Z=as.matrix(Z)
      Z=t(Z)
    }else{
      Z=as.matrix(Z)
    }
    #update membership and k
    for(i in 1 : n){
      for(j in 1 : k){
        dist[j]=wt.euc.dist.sq(X[i,],Z[j,],weight)
      }
      label[i]=which.min(dist)
      if(min(dist)>lambda_k){
        k=k+1
        Z=rbind(Z,t(as.matrix(X[i,])))
        dist=c(dist,0)
        dist1=c(dist1,0)
        label[i]=k
        for(t1 in 1:tmax1){
          for(i1 in 1 : n){
            for(j1 in 1 : k){
              dist1[j1]=wt.euc.dist.sq(X[i1,],Z[j1,],weight)
            }
            label[i1]=which.min(dist1)
          }
          for(j1 in 1:k){
            if(sum(label==j1)>0){
              I=which(label==j1)
              if(length(I)==1){
                Z[j1,]=X[I,]
              }
              else
                Z[j1,]=colMeans(X[I,])
            }
          }
        }
      }
    }
    #update centroids
    for(j in 1:k){
      if(sum(label==j)>0){
        I=which(label==j)
        if(length(I)==1){
          Z[j,]=X[I,]
        }
        else
          Z[j,]=colMeans(X[I,])
      }
    }
    #update weights
    if(t%%2==0){
      for(l in 1:p){
        D[l]=0
      }
      for(j in 1:k){
        if(sum(label==j)>0){
          I=which(label==j)
          for(i in I){
            D=D+vec.euc.dist.sq(X[i,],Z[j,],rep(1,p))
          }
        }
      }
      for(l in 1:p){
        D[l]=exp(-D[l]/lambda_w)
      }
      sum=sum(D)
      weight=D/sum
    }
    if(t>tmax){
      break
    }
    
  }
  return(list(k,weight,label))
}











