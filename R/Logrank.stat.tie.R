Logrank.stat.tie=function(x.trunc,z.trunc,d){
m=length(x.trunc)
x.grid=sort(as.numeric(levels(factor(x.trunc))))
num.x=length(x.grid)
############ truncation proportion ############
x.min=min(x.trunc)
R1=sum((x.trunc<=x.min)&(z.trunc>=x.min))
hat.c=m/R1
for(i in 2:num.x){
  v=x.grid[i]
  Nx=sum((x.trunc==v))
  Rx=sum((x.trunc<=v)&(z.trunc>=v))
  if(Rx>1){hat.c=hat.c*(1-Nx/Rx)}
}

############ Log-rank statistics ############
L0=0;L1=0;Llog=0
for(i in 1:num.x){
  xx=x.grid[i]
  z.max=max(z.trunc[as.logical((x.trunc==xx))])
  temp.z=(z.trunc>=xx)&(x.trunc<=xx)&(z.trunc<=z.max)
  z.grid=as.numeric(levels(factor(z.trunc[temp.z])))
  num.z=length(z.grid)
  for(j in 1:num.z){
    zz=z.grid[j]
    Sc=1
    z.order=sort(  as.numeric(levels(factor(z.trunc)))  )
    num.c=sum( z.order<zz )
    if(num.c>0){
      for(k in 1:num.c){
        u=z.order[k]
        Nc=sum( (1-d)[(z.trunc==u)] )
        Rc=sum((x.trunc<=u)&(z.trunc>=u))
        if(Rc>1){Sc=Sc*(1-Nc/Rc)}
      }
    }
    n11=sum((x.trunc==xx)&(z.trunc==zz)&(d==1))
    n10=sum((x.trunc==xx)&(z.trunc>=zz))
    n01=sum((x.trunc<=xx)&(z.trunc==zz)&(d==1))
    R=sum((x.trunc<=xx)&(z.trunc>=zz))
    L0=L0+n11-n10*n01/R
    hat.v=R/m/Sc
    L1=L1+hat.v*(n11-n10*n01/R)
    Llog=Llog-1/log(hat.c*hat.v)*(n11-n10*n01/R)
  }
}
c(L0,L1,Llog)
}
