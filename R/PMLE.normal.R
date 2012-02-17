PMLE.normal <-
function(l.trunc,x.trunc,testimator=FALSE){

m=length(l.trunc) ### sample size ###

############ Obtaining MLE ###############
bar.l=mean(l.trunc);bar.x=mean(x.trunc)
s.l=var(l.trunc);s.x=var(x.trunc);s.lx=cov(l.trunc,x.trunc)

l.func=function(theta){
  mul=theta[1];mux=theta[2];varl=theta[3];varx=theta[4];covlx=theta[5]
  delta=(mux-mul)/sqrt(varx+varl-2*covlx)
  prop=pnorm(delta)
  R=varl*varx-covlx^2
  D.vec=varx*(l.trunc-mul)^2-2*covlx*(l.trunc-mul)*(x.trunc-mux)+varl*(x.trunc-mux)^2
  D=mean(D.vec)/R
  -(  -m*log(prop)-m/2*log((2*pi)^2)-m/2*log(R)-m*D/2  )
}

res=nlm(l.func,p=c(bar.l,bar.x,s.l,s.x,s.lx),hessian=TRUE)
est=res$estimate
l.max=-res$minimum
Info=res$hessian/m ### Observed Fisher information ###
I_inv=solve(Info)

############ Testing quasi-independence ############
l0.func=function(theta0){
l.func(c(theta0,0))
}
res0=nlm(l0.func,p=c(bar.l,bar.x,s.l,s.x),hessian=T)
est0=res0$estimate
l0.max=-res0$minimum

LR=2*(l.max-l0.max)
est_test=c(est0,0)+(est-c(est0,0))*(LR>qchisq(0.95,df=1))

if(testimator==TRUE){est=est_test}

mul=est[1];mux=est[2];varl=est[3];varx=est[4];covlx=est[5]
se.mul=sqrt(I_inv[1,1]/m);se.mux=sqrt(I_inv[2,2]/m);se.varl=sqrt(I_inv[3,3]/m)
se.varx=sqrt(I_inv[4,4]/m);se.covlx=sqrt(I_inv[5,5]/m)

Vlx=sqrt(varx-covlx^2/varl)

###### Inclusion probability and its derivative ######
prop.func=function(mul,mux,varl,varx,covlx){
  delta=(mux-mul)/sqrt(varx+varl-2*covlx)
  pnorm(delta)
}
prop.est=prop.func(mul,mux,varl,varx,covlx)

h=0.00000001
prop_dmul=(prop.func(mul+h,mux,varl,varx,covlx)-prop.est)/h
prop_dmux=(prop.func(mul,mux+h,varl,varx,covlx)-prop.est)/h
prop_dvarl=(prop.func(mul,mux,varl+h,varx,covlx)-prop.est)/h
prop_dvarx=(prop.func(mul,mux,varl,varx+h,covlx)-prop.est)/h
prop_dcovlx=(prop.func(mul,mux,varl,varx,covlx+h)-prop.est)/h

c_dot=c(prop_dmul,prop_dmux,prop_dvarl,prop_dvarx,prop_dcovlx)

se.prop=sqrt(t(c_dot)%*%(I_inv/m)%*%c_dot)


mu_L=c(estim=mul,se=se.mul)
mu_X=c(estim=mux,se=se.mux)
var_L=c(estim=varl,se=se.varl)
var_X=c(estim=varx,se=se.varx)
cov_LX=c(estim=covlx,se=se.covlx)
prop=c(estim=prop.est,se=se.prop)
LR_test=c(LR=LR,pvalue=1-pchisq(LR,df=1))

list(mu_L=mu_L,mu_X=mu_X,var_L=var_L,var_X=var_X,cov_LX=cov_LX,
prob=prop,test=LR_test)
}
