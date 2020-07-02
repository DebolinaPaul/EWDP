#Simulation 1 and 2
data_generate=function(n,M,prob1,sigma,sigma2){
  p=dim(M)[2]
  X=matrix(0,n,p)
  k=dim(M)[1]
  label=numeric(n)
  for(i in 1:n){
    s=sample(1:k,size=1,prob=prob1)
    for(l in 1:p){
      if(M[s,l]==0){
        X[i,l]=rnorm(1,M[s,l],sigma2)
      }else{
        X[i,l]=rnorm(1,M[s,l],sigma)
      }
      
    }
    label[i]=s
  }
  ls=list(X,label)
  names(ls)=c('data','label')
  return(ls)
}

k=10
p=20
M=rand(k,p)
X=data_generate(200,M,c(0.5,0.5),0.5)
plot(X$data,col=X$label)
toss=X$label
X=X$data
plot(X,col=toss,pch=toss+13)

#heteroskedasticity
X=matrix(0,200,200)
X[1:100,1]=rnorm(100,0,0.1)
X[1:100,2]=rnorm(100,0,0.1)
X[101:200,1]=rnorm(100,3.2,1)
X[101:200,2]=rnorm(100,3.2,1)
#X[201:300,1]=rnorm(100,10,2)
#X[201:300,2]=rnorm(100,-10,2)
for(l in 3:200){
  X[,l]=rnorm(200,0,2)
}
toss=c(rep(1,100),rep(2,100))
plot(X,col=toss)
write.csv(X,file = 'cor.csv')

#Correlated
X=matrix(0,200,200)
X[1:100,1:2]=mvtnorm::rmvnorm(100,c(0,0),matrix(c(1,0.5,0.5,1),nrow = 2))
X[101:200,1:2]=mvtnorm::rmvnorm(100,c(7,7),matrix(c(1,-0.5,-0.5,1),nrow = 2))
#X[,3:4]=mvtnorm::rmvnorm(200,c(0,0),matrix(c(2,0.5,0.5,2),nrow = 2))
for(l in 3:200){
  X[,l]=rnorm(200,0,5)
}
toss=c(rep(1,100),rep(2,100))
plot(X,col=toss,pch=19)

