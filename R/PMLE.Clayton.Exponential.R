PMLE.Clayton.Exponential <-
function(l.trunc,x.trunc,GOF=TRUE,Err=3,alpha_max=20,alpha_min=10^-4){
l=l.trunc
x=x.trunc


A_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  FL = 1-exp(-lamdaL*l)
  
  FL^(-alpha)+FX^(-alpha)-1
}

A_a_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  FL = 1-exp(-lamdaL*l)
  
  -log(FL)/(FL^alpha)-log(FX)/(FX^alpha)
}

A_L_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FL = 1-exp(-lamdaL*l)
  
  -alpha*l*exp(-lamdaL*l)/(FL^(alpha+1))
}

A_X_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  
  -alpha*x*exp(-lamdaX*x)/(FX^(alpha+1))
}

A_aa_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  FL = 1-exp(-lamdaL*l)
  
  (log(FL))^2/(FL^alpha)+(log(FX))^2/(FX^alpha)
}

A_LL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FL = 1-exp(-lamdaL*l)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  ((-alpha-1)*l*exp(-lamdaL*l)/FL-l)*A_L
}

A_XX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  ((-alpha-1)*x*exp(-lamdaX*x)/FX-x)*A_X
}

A_aL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FL = 1-exp(-lamdaL*l)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  (1/alpha-log(FL))*A_L
}

A_aX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  (1/alpha-log(FX))*A_X
}

######################################################################################################################################

B_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U=1-(1-u)^(lamdaL/lamdaX) 
  
  U^(-alpha)+u^(-alpha)-1
}

B_a_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  
  -log(U)/U^alpha-log(u)/u^alpha
}

B_L_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  
  alpha*(1-U)*log(1-u)/(lamdaX*U^(alpha+1))
}

B_X_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  
  -alpha*lamdaL*(1-U)*log(1-u)/(lamdaX^2*U^(alpha+1))
}

B_aa_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  
  (-log(U))^2/(U^alpha)+(-log(u))^2/u^alpha
}

B_LL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  B_L = B_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  (((alpha+1)*(1-U))*log(1-u)/(lamdaX*U)+log(1-u)/lamdaX)*B_L
}

B_XX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  B_X = B_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  ((alpha+1)*(1-U)*log(1-u)/U+log(1-u)+2*lamdaX/lamdaL)*(-lamdaL/lamdaX^2)*B_X
}

B_aL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_L = B_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  (1/alpha-log(U))*B_L
}

B_aX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_X = B_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  (1/alpha-log(U))*B_X
}

B_LX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  U = 1-(1-u)^(lamdaL/lamdaX)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_X = B_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  ((alpha+1)*(1-U)*log(1-u)/(lamdaX*U)+1/lamdaL+log(1-u)/lamdaX)*B_X
}

########################################################################################################################################

###truncated part
H_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  u^(-alpha-1)*B^(-1/alpha-1)
}

###first_order
H_a_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_a = B_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  (u^(-alpha-1)*B^(-1/alpha-1)*(-log(u)+log(B)/alpha^2+(-1/alpha-1)*B_a/B))*alpha
}

H_L_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_L = B_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  ((-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_L)*lamdaL
}

H_X_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_X = B_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  
  ((-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_X)*lamdaX
}

###second
H_aa_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_a = B_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_aa = B_aa_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  H_a = H_a_func(u)
  
  H_a+((-log(u)+log(B)/alpha^2+(-1/alpha-1)*B_a/B)*H_a/alpha+
         u^(-alpha-1)*B^(-1/alpha-1)*(2*B_a/(alpha^2*B)-2*log(B)/alpha^3+
                                        (-1/alpha-1)*(B_aa*B-B_a^2)/B^2))*alpha^2
}

H_aL_func= function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_a = B_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_L = B_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_aL = B_aL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  H_a = H_a_func(u)
  
  ((-1/alpha-1)*B_L/B*H_a/alpha+
    u^(-alpha-1)*B^(-1/alpha-1)*(B_L/(B*alpha^2)+
                                   (-1/alpha-1)*(B_aL*B-B_L*B_a)/B^2))*lamdaL*alpha
}

H_aX_func= function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_a = B_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_X = B_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_aX = B_aX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  H_a = H_a_func(u)
  
  ((-1/alpha-1)*B_X/B*H_a/alpha+
    u^(-alpha-1)*B^(-1/alpha-1)*(B_X/(B*alpha^2)+
                                   (-1/alpha-1)*(B_aX*B-B_X*B_a)/B^2))*lamdaX*alpha
}

H_LL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_L = B_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_LL = B_LL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  H_L = H_L_func(u)
  
  H_L+((-1/alpha-2)*B_L/B*H_L/lamdaL+(-1/alpha-1)*B_LL/(u^(alpha+1)*B^(1/alpha+2)))*lamdaL^2
}

H_XX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_X = B_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_XX = B_XX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  H_X = H_X_func(u)
  
  H_X+((-1/alpha-2)*B_X/B*H_X/lamdaX+(-1/alpha-1)*B_XX/(u^(alpha+1)*B^(1/alpha+2)))*lamdaX^2
}

H_LX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  B = B_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_L = B_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_X = B_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  B_LX = B_LX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,u)
  H_X = H_X_func(u)
  
  ((-1/alpha-2)*B_L/B*H_X/lamdaX+(-1/alpha-1)*B_LX/(u^(alpha+1)*B^(1/alpha+2)))*lamdaL*lamdaX
}

###complete part
clayton_complete = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  FL = 1-exp(-lamdaL*l)
  fX = lamdaX*exp(-lamdaX*x)
  fL = lamdaL*exp(-lamdaL*l)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  fL*fX*(1+alpha)*(FL*FX)^(-alpha-1)/A^(1/alpha+2)
}

pdf_a_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  FL = 1-exp(-lamdaL*l)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_a = A_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  (1/(1+alpha)-log(FL)-log(FX)+1/alpha^2*log(A)+(-1/alpha-2)*A_a/A)*alpha
}

alpha_tuda = 1
lamdaL_tuda = 4
lamdaX_tuda = 5
pdf_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,2,1)
a=log(clayton_complete(alpha_tuda,lamdaL_tuda,lamdaX_tuda,2,1))
alpha_tuda = alpha_tuda+10^-8
b=log(clayton_complete(alpha_tuda,lamdaL_tuda,lamdaX_tuda,2,1))
(b-a)/10^-8

pdf_L_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FL = 1-exp(-lamdaL*l)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  ((-alpha-1)*l*exp(-lamdaL*l)/FL+(-1/alpha-2)*A_L/A+1/lamdaL-l)*lamdaL
}

pdf_X_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  ((-alpha-1)*x*exp(-lamdaX*x)/FX+(-1/alpha-2)*A_X/A+1/lamdaX-x)*lamdaX
}

##########
pdf_aa_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_a = A_a_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_aa = A_aa_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x) 
  pdf_a = pdf_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  pdf_a+(-1/(1+alpha)^2-2/alpha^3*log(A)+2/alpha^2*A_a/A+
           (-1/alpha-2)*(A_aa*A-A_a^2)/A^2)*alpha^2
}

pdf_LL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FL = 1-exp(-lamdaL*l)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_L = A_L_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_LL = A_LL_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x) 
  pdf_L = pdf_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  pdf_L+((alpha+1)*(l^2*exp(-lamdaL*l))/FL^2+(-1/alpha-2)*(A_LL*A-A_L^2)/A^2
         -1/lamdaL^2)*lamdaL^2
}

pdf_XX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_X = A_X_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_XX = A_XX_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x) 
  pdf_X = pdf_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  pdf_X+((alpha+1)*(x^2*exp(-lamdaX*x))/FX^2+(-1/alpha-2)*(A_XX*A-A_X^2)/A^2-
           1/lamdaX^2)*lamdaX^2
}

pdf_aL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FL = 1-exp(-lamdaL*l)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_L = A_L_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_a = A_a_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_aL = A_aL_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x) 
  
  (-l*exp(-lamdaL*l)/FL+A_L/(alpha^2*A)+(-1/alpha-2)*(A_aL*A-A_a*A_L)/A^2)*
    alpha*lamdaL
}

pdf_aX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  FX = 1-exp(-lamdaX*x)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_a = A_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_aX = A_aX_func (alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x) 
  
  (-x*exp(-lamdaX*x)/FX+A_X/(alpha^2*A)+(-1/alpha-2)*(A_aX*A-A_a*A_X)/A^2)*
    alpha*lamdaX
}

pdf_LX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  ((-1/alpha-2)*(-A_X*A_L)/A^2)*lamdaL*lamdaX
}

hessian_func = function(er){    
  prob = integrate(H_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_a = integrate(H_a_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_L = integrate(H_L_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_X = integrate(H_X_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_aa = integrate(H_aa_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_LL = integrate(H_LL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_XX = integrate(H_XX_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_aL = integrate(H_aL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_aX = integrate(H_aX_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_LX = integrate(H_LX_func, lower = 0, upper = 1,rel.tol = er)$value
  
  pdf_a = pdf_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_L = pdf_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_X = pdf_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_aa = pdf_aa_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_LL = pdf_LL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_XX = pdf_XX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_aL = pdf_aL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_aX = pdf_aX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_LX = pdf_LX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  hessian = matrix(c(
    -n/prob^2*(prob_aa*prob-prob_a^2)+sum(pdf_aa),
    -n/prob^2*(prob_aL*prob-prob_a*prob_L)+sum(pdf_aL),
    -n/prob^2*(prob_aX*prob-prob_a*prob_X)+sum(pdf_aX),
    -n/prob^2*(prob_aL*prob-prob_a*prob_L)+sum(pdf_aL),
    -n/prob^2*(prob_LL*prob-prob_L^2)+sum(pdf_LL),
    -n/prob^2*(prob_LX*prob-prob_X*prob_L)+sum(pdf_LX),
    -n/prob^2*(prob_aX*prob-prob_a*prob_X)+sum(pdf_aX),
    -n/prob^2*(prob_LX*prob-prob_X*prob_L)+sum(pdf_LX),
    -n/prob^2*(prob_XX*prob-prob_X^2)+sum(pdf_XX)
  ),3,3,byrow = TRUE)
  hessian
}

score_func = function(er){    
  prob = integrate(H_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_a = integrate(H_a_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_L = integrate(H_L_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_X = integrate(H_X_func, lower = 0, upper = 1,rel.tol = er)$value
  
  
  pdf_a = pdf_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_L = pdf_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  pdf_X = pdf_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,l,x)
  
  score = c(-n*prob_a/prob+sum(pdf_a),-n*prob_L/prob+sum(pdf_L),-n*prob_X/prob+sum(pdf_X))
  
  
  score
}

rel_tol_func = function(CI){
  
  e =.Machine$double.eps^0.25*10^sample(-4:2,1) 
  repeat{
    EE = ( 
      try(integrate(H_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_a_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_L_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_aa_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_LL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_aL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_LX_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 )
    if( EE == T | CI==3 ){break}
    if( EE == F ){ e =.Machine$double.eps^0.25*10^CI ; CI=CI+1} 
    
  }
  e
}

O_pdf_func = function(alpha,lamdaL,lamdaX,l,x){
  (1+alpha)*lamdaL*exp(-lamdaL*l)*lamdaX*exp(-lamdaX*x)*((1-exp(-lamdaL*l))*
                                                           (1-exp(-lamdaX*x)))^(-alpha-1)*((1-exp(-lamdaL*l))^-alpha+
                                                                                             (1-exp(-lamdaX*x))^-alpha-1)^(-1/alpha-2)
}
O_h_func = function(u){
  u^(-alpha-1)*((1-(1-u)^(lamdaL/lamdaX))^-alpha+u^-alpha-1)^(-1-1/alpha)
}

l=l.trunc
x=x.trunc
n=length(x)
par_A = c()
par_L = c()
par_X = c()
par_new = matrix(c(NA,NA,NA),3,1)
par_old = matrix(c(NA,NA,NA),3,1)

tau=cor(l,x,method = "kendall")
par_old[1,1] = log(2*tau/(1-tau))
par_old[2,1] = log(1/mean(l))
par_old[3,1] = log(1/mean(x))

alpha_tuda = par_old[1,1];lamdaL_tuda = par_old[2,1];lamdaX_tuda = par_old[3,1]
par_A[1] = par_old[1,1];par_L[1] = par_old[2,1];par_X[1] = par_old[3,1]
AI = 1
BI = 1
e =rel_tol_func(-5)
repeat{
  par_new = par_old-solve(hessian_func(e))%*%score_func(e)
  alpha_tuda = par_new[1,1]
  lamdaL_tuda = par_new[2,1]
  lamdaX_tuda = par_new[3,1]
  error1 = abs(exp(par_new[1,1])-exp(par_old[1,1]))
  error2 = abs(exp(par_new[2,1])-exp(par_old[2,1]))
  error3 = abs(exp(par_new[3,1])-exp(par_old[3,1]))    
  Error1 = max(error1,error2,error3)
  Error2 = ( Error1 > Err | exp(par_new[1,1]) < alpha_min | 
               exp(par_new[2,1]) < 10^-10 | exp(par_new[3,1]) < 10^-10 |
               exp(par_new[1,1]) > alpha_max | AI/50 == floor(AI/50))    
  e = rel_tol_func(-5)
  Error3 = try(max(eigen(hessian_func(e))$value),silent = T)
  if( Error1 <10^-4 & Error3 < 0 ) {break}
  if( Error2 == T ){
    repeat{
      par_new[1,1] = log(2*tau/(1-tau)*exp(runif(1,-1.5,1.5)))
      par_new[2,1] = log(1/mean(l)*exp(runif(1,-0.5,0.5)))
      par_new[3,1] = log(1/mean(x)*exp(runif(1,-0.5,0.5)))    
      BI=BI+1
      if( exp(par_new[1,1])<alpha_max & exp(par_new[1,1])>alpha_min ){break}
    }
  }
  par_old = par_new
  AI = AI+1 
  alpha_tuda = par_old[1,1]
  lamdaL_tuda = par_old[2,1]
  lamdaX_tuda = par_old[3,1]
  e = rel_tol_func(-5)
  par_A[AI] = par_old[1,1];par_L[AI] = par_old[2,1];par_X[AI] = par_old[3,1]
}
alpha_tuda = par_new[1,1]
lamdaL_tuda = par_new[2,1]
lamdaX_tuda = par_new[3,1]

alpha = exp( par_new[1,1])
lamdaL = exp( par_new[2,1])
lamdaX = exp( par_new[3,1])
mean_X = 1/lamdaX

TT = solve(-hessian_func(e))
se_alpha=exp(alpha_tuda)*sqrt(TT[1,1])
se_lamdaL=exp(lamdaL_tuda)*sqrt(TT[2,2])
se_lamdaX=exp(lamdaX_tuda)*sqrt(TT[3,3])
se_mu=sqrt(TT[3,3])/lamdaX
alpha_res = c(Estimate=alpha,SE=se_alpha)
lambda_L_res = c(Estimate=lamdaL,SE=se_lamdaL)
lambda_X_res = c(Estimate=lamdaX,SE=se_lamdaX)
mu_res = c(Estimate=mean_X,SE=se_mu)
LL=-n*log(integrate(O_h_func, lower = 0, upper = 1)$value)+
  sum(log(O_pdf_func(alpha,lamdaL,lamdaX,l,x)))
AIC_res=-2*LL+2*5
BIC_res=-2*LL+5*log(n)

C.test=K.test=NULL
F_par=F_emp=prop=NULL

if(GOF==TRUE){
nuL=1;nuX=1
  
F.func=function(ll,xx){
  Fll=1-exp(-lamdaL*ll^nuL)
  Fxl=1-exp(-lamdaX*ll^nuX)
  Fxx=1-exp(-lamdaX*xx^nuX)
  F1 = integrate(H_func, lower = 0, upper = Fxl)$value
  HH_func = function(u){
    B=Fll^(-alpha)+u^(-alpha)-1
    u^(-alpha-1)*B^(-1/alpha-1)
  }
  F2 = integrate(HH_func, lower = Fxl, upper = Fxx)$value
  F1+F2
}

prop=F.func(Inf,Inf)
F_par=F_emp=numeric(n)
for(i in 1:n){
  F_par[i]=F.func(l[i],x[i])/prop
  F_emp[i]=mean( (l<=l[i])&(x<=x[i]) )
}
C.test=sum( (F_emp-F_par)^2 )
K.test=max( abs( F_emp-F_par ) )

plot(F_emp,F_par,xlab="F_empirical",ylab="F_parametric",xlim=c(0,1),ylim=c(0,1))
lines(x = c(0,1), y = c(0,1))

}


list(n=n,alpha = alpha_res,lambda_L = lambda_L_res,lambda_X = lambda_X_res
     ,mean_X = mu_res,logL=round(LL,2),AIC=AIC_res,BIC=BIC_res,Iteration=AI,
     c=prop,C=C.test,K=K.test,F_empirical=F_emp,F_parametric=F_par)
}

