Survival estimation Coxme
================
Santiago Bohórquez
22 de junio de 2016

Estimation of baseline survival
-------------------------------

Using Efron for tied deaths.

``` r
require(coxme)
```

    ## Loading required package: coxme

    ## Loading required package: survival

    ## Loading required package: bdsmatrix

    ## 
    ## Attaching package: 'bdsmatrix'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

``` r
efron_me<-function(alpha,cox_me,i){
  
  if (ncol(cox_me$y)==3) {R_m<-which(cox_me$y[,1]<i &
                                    cox_me$y[,2]>=i)}
  else {R_m<-which(cox_me$y[,1]>=i)}
  D_j<-intersect(R_m, which(cox_me$y[,ncol(cox_me$y)]==1 &
                              cox_me$y[,(ncol(cox_me$y)-1)]==i ))
  
  rhs<-sum(exp(cox_me$linear.predictor[R_m]))

  lhs_exp<-exp(cox_me$linear.predictor[D_j])
  lhs<-sum(lhs_exp/(1-alpha^lhs_exp))
  
  return(abs(rhs-lhs))
}

alpha_j<-function(cox_me,j){
  op<-optim(0.5,efron_me, cox_me=cox_me, i=j,  method = "Brent",
            lower=0, upper=1, control=list(abstol=1e-8))
  #op<-nleqslv(0.5, efron_me,cox_me=cox31,i=1) #sligthly slower
  return(op$par)
}

S0<-function(cox_me){
  t<-max(cox_me$y[,(ncol(cox_me$y)-1)])
  alpha<-sapply(1:t,alpha_j,cox_me=cox_me)
  return(cumprod(alpha))
}

survfit.coxme<-function(cox_me,newdata=NULL){
  survfit<-list()
  survfit$n<-cox_me$n[2]
  if (length(newdata)==0) {
    S0_est<-S0(cox_me)
    survfit$surv<-S0_est^(exp(sum(cox_me$means*cox_me$coefficients))) 
    return(survfit)
  }
  else {
    S0_est<-S0(cox_me)
    data_est<-data.matrix(newdata[,names(cox_me$coefficients)])%*%cox_me$coefficients
    survfit$surv<-sapply(data_est,function(x,y) y^(exp(x)),y=S0_est)
    return(survfit)
  }
}
# 
# print.survfit.coxme<-function(x,...){
#   cat("Surv:\n")
#   print(x$surv)
# }
```

``` r
attach(lung)

lung_coxme<-coxme(Surv(time,status)~age+sex+(1|inst))
# newdata<-cbind.data.frame(eneindex=c(1,2),unemp=c(500,500))
summary(lung_coxme)
```

    ## Cox mixed-effects model fit by maximum likelihood
    ## 
    ##   events, n = 164, 227 (1 observation deleted due to missingness)
    ##   Iterations= 5 22 
    ##                     NULL Integrated    Fitted
    ## Log-likelihood -744.7999  -737.8121 -737.8009
    ## 
    ##                   Chisq   df          p  AIC   BIC
    ## Integrated loglik 13.98 3.00 0.00293850 7.98 -1.32
    ##  Penalized loglik 14.00 2.01 0.00092769 9.97  3.74
    ## 
    ## Model:  Surv(time, status) ~ age + sex + (1 | inst) 
    ## Fixed coefficients
    ##            coef exp(coef)    se(coef)     z      p
    ## age  0.01703526 1.0171812 0.009233008  1.85 0.0650
    ## sex -0.51167566 0.5994902 0.167681629 -3.05 0.0023
    ## 
    ## Random effects
    ##  Group Variable  Std Dev      Variance    
    ##  inst  Intercept 9.110830e-03 8.300722e-05

``` r
jaja<-survfit(lung_coxme)

jum<-survfit(coxph(Surv(time,status)~age+sex))
```

Probblems so far:

-   Factors are not supported
-   Calculates surv at all times, better to only calculate given a list or where there are events.
-   No confidence intervals
-   Only returns n and surv
