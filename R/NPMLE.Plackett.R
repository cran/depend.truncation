NPMLE.Plackett <-
function(x.trunc,y.trunc,x.fix=median(x.trunc),y.fix=median(y.trunc),
Fx.plot=TRUE){

  m=length(x.trunc)
  x.ox=sort(x.trunc);y.oy=sort(y.trunc)
  dHx1=1;Hxn=0;Ay_1=0;dAyn=1 ### initial constraint ####

  l.func=function(dL){
    dHx=c(dHx1,exp(dL[1:(m-1)]));dAy=c(exp(dL[m:(2*m-2)]),dAyn)
    Hx=c(rev(cumsum(rev(dHx)))[-1],Hxn);Ay_=c(Ay_1,cumsum(dAy)[-m])
    alpha=exp(dL[2*m-1])  

    prop=0
    for(i in 1:m){
      temp_y=(y.oy>=x.ox[i])
      dHxi=dHx[i];dAyi=dAy[temp_y]
      Fxi=exp(-Hx[i]);Sy_i=exp(-Ay_[temp_y])
      B=1+(alpha-1)*(Fxi+Sy_i)
      C=(alpha-1)*Fxi*Sy_i
      c11=alpha*(B-2*C)/(B^2-4*alpha*C)^(3/2)
      prop=prop+sum(c11*Sy_i*dAyi*Fxi*dHxi)
    }
    l=-m*log(prop)
    for(i in 1:m){
      x_num=sum(x.ox<=x.trunc[i]);y_num=sum(y.oy<=y.trunc[i])
      Hxi=Hx[x_num];Ay_i=Ay_[y_num]
      dHxi=dHx[x_num];dAyi=dAy[y_num]
      Fxi=exp(-Hxi);Sy_i=exp(-Ay_i)
      B=1+(alpha-1)*(Fxi+Sy_i)
      C=(alpha-1)*Fxi*Sy_i
      l=l+log(alpha)+log(B-2*C)-3/2*log(B^2-4*alpha*C)-Hxi-Ay_i+log(dHxi)+log(dAyi)
    }
    -l
  }

  ##### Initial value by the product-limit estimates ######
  dHx_lyn=dAy_lyn=numeric(m)
  for(i in 1:m){
    dHx_lyn[i]=1/sum(  (x.trunc<=x.ox[i])& (y.trunc>=x.ox[i])  )
    dAy_lyn[i]=1/sum(  (x.trunc<=y.oy[i])& (y.trunc>=y.oy[i])  )
  }
  dHx_lyn=dHx_lyn[-1]
  dAy_lyn=dAy_lyn[-m]
  dL_lyn=c(log(c(dHx_lyn,dAy_lyn)),log(1.0001))

  res=nlm(l.func,p=dL_lyn,hessian=TRUE)
  dL=res$estimate;conv=res$code
  alpha=exp(dL[2*m-1])
  dHx=c(dHx1,exp(dL[1:(m-1)]))
  dAy=c(exp(dL[m:(2*m-2)]),dAyn)
  Hx_=rev(cumsum(rev(dHx)))
  Ay=cumsum(dAy)
  Fx=exp(-Hx_)
  Sy=exp(-Ay)
  Hx=c(Hx_[-1],Hxn)
  Ay_=c(Ay_1,Ay[-m])
  V_dL=solve(res$hessian)
  V=diag(exp(dL))%*%V_dL%*%diag(exp(dL))
  SE_alpha=sqrt( V_dL[2*m-1,2*m-1] )

  Hx_fix=Ay_fix=numeric(length(x.fix))
  SE_Hx=SE_Ay=numeric(length(x.fix))

  for(i in 1:length(x.fix)){
    temp.x=c(x.ox>=x.fix[i])
    Hx_fix[i]=Hx_[max(sum(x.ox<=x.fix[i]),1)]
    SE_Hx[i]=sqrt( (temp.x[-1])%*%V[1:(m-1),1:(m-1)]%*%(temp.x[-1]) )
  }   

  for(i in 1:length(y.fix)){
    temp.y=c(y.oy<=y.fix[i])
    Ay_fix[i]=Ay[max(sum(y.oy<=y.fix[i]),1)]
    SE_Ay[i]=sqrt( (temp.y[-m])%*%V[m:(2*m-2),m:(2*m-2)]%*%(temp.y[-m]) )
  }   

  Fx_fix=exp(-Hx_fix);Sy_fix=exp(-Ay_fix)
  SE_Fx=Fx_fix*SE_Hx;SE_Sy=Sy_fix*SE_Ay
  if(Fx.plot==TRUE){plot(x.ox,Fx,type="s",xlab="x",ylab="Fx(x)")}
  list(alpha=alpha,alpha_se=SE_alpha,Hx=Hx_fix,Hx_se=SE_Hx,Ay=Ay_fix,Ay_se=SE_Ay,
Fx=Fx_fix,Fx_se=SE_Fx,Sy=Sy_fix,Sy_se=SE_Sy,convergence=conv)
}
