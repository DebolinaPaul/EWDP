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

data(iris)
X=iris
X=data.matrix(X)
toss=X[,5]
X=X[,-5]
l=ERAkmeans(X,1,0.555,50)
l[[1]]
plot(X,col=l[[3]])
compare(l[[3]],toss,'nmi')
compare(l[[3]],toss,'adjusted.rand')

#kmeans
l=kmeans(X,28)
compare(l[[1]],toss,'nmi')
compare(l[[1]],toss,'adjusted.rand')
plot(X,col=l[[1]],pch=19)

setwd("E:/Educational/r works/k.means/Datasets/classification keel")
X=read.csv("abalone.csv",header = F)
X=data.matrix(X)
p=dim(X)[2]
toss=X[,p]
X=X[,-p]
table(toss)
#Scale
p=dim(X)[2]
for(j in 1:p){
  X[,j]=(X[,j]-sum(X[,j]))/sd(X[,j])
}
plot(X,col=toss+1,pch=toss+13)
#Run
l=ERAkmeans(X,0.01,0.02,50)
length(unique(l[[3]]))
l[[1]]
plot(X,col=l[[3]],pch=19)
compare(l[[3]],toss,'nmi')
compare(l[[3]],toss,'adjusted.rand')
l[[2]]
hist(l[[2]])
(ind=which(l[[2]]>0.002))
plot((1:p),l[[2]])

X=readMat('warpAR10p.mat')
toss=X$Y
X=X$X
Y=tsne::tsne(dist(X))
plot(Y,col=l[[3]],pch=l[[3]]+15)

#tsne
Z=tsne::tsne(dist(X))
plot(Z,col=toss,pch=19)
plot(Z,col=l[[3]],pch=19)

#ggplot2
Y=cbind(X,l1[[3]])
Y=data.frame(Y)
names(Y)
Y$X11=as.factor(Y$X11)
g=ggplot(data=Y,aes(x=X1,y=X2,col=X11))+geom_point(size=2)
g=g+labs(x='Dimension 1',y='Dimension 2')
g=g+scale_colour_brewer(palette='Set1')
g=g+scale_colour_wsj('colors6')
g=g+scale_f_brewer(palette="Set1")
g=g+theme(
  # axis.title.x=element_text(size=rel(5)),
  # axis.title.y=element_text(size=rel(5)),
  axis.text = element_text(size=12),
  axis.title.x = element_text( size=20),
  axis.title.y = element_text( size=20),
  legend.position = "none"
)
g

pdf("d1_10_dpmeans.pdf",height=4,width=4)
g
dev.off()

#Microarray
setwd("E:\\Educational\\r works\\k.means\\Datasets\\microarray")
X=read.table("leukemia.x.txt")
X=data.matrix(X)
X=t(X)
toss=read.table("leukemia.y.txt")
toss=data.matrix(toss)
table(toss)

#weight analysis

wdp=matrix(0,20,20)
wk=matrix(0,20,20)
wg=matrix(0,20,20)
for(m in 1:20){
  X=data_generate(100,M,rep(1/k,k),0.015,1)
  X=X$data
  l1=ERAkmeans(X,1,1,50)
  wdp[,10]=l1[[2]]
  l2=wkmeans(X,10,2,10)
  wk[,10]=sample(20,10)
  wk[,10]=wk[,10]/sum(wk[,10])
  l3=wgmeans(X,2,0.0005,30)
  m=4
  wg[,m]=sample(20,10)
  wg[,m]=wg[,m]/sum(wg[,m])
}
pdf("w_wg.pdf",height=4,width=4)
boxplot(t(wg[,1:4]),col='orange')
dev.off()










