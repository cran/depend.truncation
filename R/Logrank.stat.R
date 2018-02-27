Logrank.stat=function(x.trunc,z.trunc,d){
m=length(x.trunc)
#######Lynden-Bell's estimator######
t=c(x.trunc,z.trunc) ;t.o=t[order(t)] ;d.o=c(rep(0,m),rep(1,m))[order(t)]
dd.o=c(rep(0,m),d)[order(t)] ;r.diag=numeric(2*m)
for(i in 1:(2*m)){
r.diag[i]=sum( (x.trunc<=t.o[i])&(z.trunc>=t.o[i]) )
}
sc.diag=cumprod( 1-d.o*(1-dd.o)/r.diag*(r.diag>1) )
sc=(sc.diag[d.o==1])[rank(z.trunc)]

#######estimating truncation proportion######
x.o=sort(x.trunc)[-1]
r.diag=numeric(m-1)
for(i in 1:(m-1)){
r.diag[i]=sum((x.trunc<=x.o[i])&(z.trunc>=x.o[i]))
}
hat.c=m*prod((1-1/r.diag)[r.diag>1])

######summarizing data(making risk set on grids)#####
r.point=sckm.point=r.grid=sckm.grid=0
for(i in 1:m){
x=x.trunc[i] ;z=z.trunc[i]
if(d[i]==1){
r.point=c(r.point,sum((x.trunc<=x)&(z.trunc>=z)))
sckm.point=c(sckm.point,sc[i])
}
index=(x.trunc<=x)&(z.trunc>x)&(z.trunc<=z)&(d==1)
num=sum(index)
if(num>0){
  r.vec=numeric(num); sckm.vec<-numeric(num)
  for(j in 1:num){
  zz=z.trunc[index][j]
  r.vec[j]=sum( (x.trunc<=x)&(z.trunc>=zz) )
  sckm.vec[j]=(sc[index])[j]
  }
  r.grid=c(r.grid,r.vec); sckm.grid=c(sckm.grid,sckm.vec)
}
}
r.point=r.point[-1] ;sckm.point=sckm.point[-1] ;r.grid=r.grid[-1]
sckm.grid=sckm.grid[-1]

######computing the logrank statistics#######
L0=length(r.point)-sum(1/r.grid) # Logrank stat. rho=0
L1=sum(r.point/sckm.point)/m-sum(1/sckm.grid)/m # Gehan's stat. rho=1
Llog=sum(1/log(hat.c*r.point/sckm.point/m))-sum(1/log(hat.c*r.grid/sckm.grid/m)/r.grid)
c(L0,L1,Llog)
}
