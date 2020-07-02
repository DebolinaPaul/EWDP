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



