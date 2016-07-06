Survival estimation Coxme
================
Santiago Bohórquez
July 6, 2016

Estimation of baseline survival
-------------------------------

Using Efron for tied deaths.

``` r
rm(list=ls(all=TRUE))

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

S0<-function(cox_me,t){
  alpha<-sapply(t,alpha_j,cox_me=cox_me)
  return(cumprod(alpha))
}

survfit.coxme<-function(cox_me,newdata=NULL){
  survfit<-list()
  survfit$n<-cox_me$n[2]
  events<-table(cox_me$y[,(ncol(cox_me$y)-1)],cox_me$y[,ncol(cox_me$y)])
  survfit$time<-as.numeric(rownames(events))
  rownames(events)<-c()
  survfit$n.risk<-cox_me$n[2]-cumsum(rowSums(events))
  survfit$n.event<-events[,2]
  survfit$n.censor<-events[,1]
  S0_est<-S0(cox_me,survfit$time)
  if (length(newdata)==0) {
    survfit$surv<-S0_est^(exp(sum(cox_me$means*cox_me$coefficients))) 
  }
   
  else {
    data_est<-data.matrix(newdata[,names(cox_me$coefficients)])%*%cox_me$coefficients
    survfit$surv<-sapply(data_est,function(x,y) y^(exp(x)),y=S0_est)
    survfit$time<-unique(cox_me$y[which(cox_me$y[,ncol(cox_me$y)]==1),
                     (ncol(cox_me$y)-1)])
  }
  survfit$type<-attr(cox_me$y,"type")
  survfit$cumhaz<- -log(survfit$surv)
  survfit$call<-match.call()
  attr(survfit,"class")<-c("survfit.cox","survfit")
 return(survfit)
  
}
```

``` r
attach(lung)

lung_coxme<-coxme(Surv(time,status)~age+sex+(1|inst))

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
summary(jaja)
```

    ## Call: survfit.coxme(cox_me = lung_coxme)
    ## 
    ##  time n.risk n.event survival
    ##     5    226       1   0.9958
    ##    11    223       3   0.9831
    ##    12    222       1   0.9789
    ##    13    220       2   0.9704
    ##    15    219       1   0.9662
    ##    26    218       1   0.9619
    ##    30    217       1   0.9576
    ##    31    216       1   0.9533
    ##    53    214       2   0.9447
    ##    54    213       1   0.9404
    ##    59    212       1   0.9361
    ##    60    210       2   0.9275
    ##    61    209       1   0.9232
    ##    62    208       1   0.9188
    ##    65    206       2   0.9102
    ##    71    205       1   0.9059
    ##    79    204       1   0.9016
    ##    81    202       2   0.8929
    ##    88    200       2   0.8843
    ##    92    198       1   0.8800
    ##    93    197       1   0.8756
    ##    95    195       2   0.8670
    ##   105    193       1   0.8626
    ##   107    191       2   0.8539
    ##   110    190       1   0.8495
    ##   116    189       1   0.8451
    ##   118    188       1   0.8407
    ##   122    187       1   0.8363
    ##   131    186       1   0.8319
    ##   132    184       2   0.8231
    ##   135    183       1   0.8188
    ##   142    182       1   0.8144
    ##   144    181       1   0.8100
    ##   145    179       2   0.8012
    ##   147    178       1   0.7968
    ##   153    177       1   0.7924
    ##   156    175       2   0.7836
    ##   163    172       3   0.7703
    ##   166    170       2   0.7615
    ##   167    169       1   0.7571
    ##   170    168       1   0.7527
    ##   175    164       1   0.7482
    ##   176    163       1   0.7437
    ##   177    161       1   0.7392
    ##   179    159       2   0.7301
    ##   180    158       1   0.7255
    ##   181    156       2   0.7164
    ##   182    155       1   0.7119
    ##   183    154       1   0.7073
    ##   186    152       1   0.7027
    ##   189    150       1   0.6981
    ##   194    147       1   0.6934
    ##   197    144       1   0.6887
    ##   199    143       1   0.6839
    ##   201    141       2   0.6744
    ##   202    139       1   0.6697
    ##   207    137       1   0.6649
    ##   208    136       1   0.6601
    ##   210    135       1   0.6553
    ##   212    133       1   0.6505
    ##   218    132       1   0.6457
    ##   222    129       1   0.6408
    ##   223    128       1   0.6358
    ##   226    124       1   0.6308
    ##   229    123       1   0.6257
    ##   230    122       1   0.6206
    ##   239    118       2   0.6103
    ##   245    115       1   0.6051
    ##   246    114       1   0.5998
    ##   267    110       1   0.5945
    ##   268    109       1   0.5891
    ##   269    107       1   0.5838
    ##   270    106       1   0.5784
    ##   283    102       1   0.5728
    ##   284    100       1   0.5672
    ##   285     98       2   0.5558
    ##   286     97       1   0.5501
    ##   288     96       1   0.5445
    ##   291     95       1   0.5387
    ##   293     92       1   0.5329
    ##   301     88       1   0.5270
    ##   303     86       1   0.5209
    ##   305     85       1   0.5148
    ##   306     84       1   0.5087
    ##   310     82       2   0.4963
    ##   320     80       1   0.4901
    ##   337     78       1   0.4838
    ##   340     77       1   0.4776
    ##   345     76       1   0.4713
    ##   348     75       1   0.4651
    ##   350     74       1   0.4589
    ##   351     73       1   0.4528
    ##   353     71       2   0.4404
    ##   361     69       1   0.4341
    ##   363     67       2   0.4216
    ##   364     65       1   0.4153
    ##   371     63       2   0.4026
    ##   387     59       1   0.3961
    ##   390     58       1   0.3895
    ##   394     57       1   0.3829
    ##   426     54       1   0.3760
    ##   428     53       1   0.3691
    ##   429     52       1   0.3621
    ##   433     51       1   0.3551
    ##   442     50       1   0.3482
    ##   444     48       1   0.3412
    ##   450     47       1   0.3340
    ##   455     46       1   0.3269
    ##   457     45       1   0.3197
    ##   460     43       1   0.3123
    ##   473     42       1   0.3048
    ##   477     41       1   0.2974
    ##   519     38       1   0.2896
    ##   520     37       1   0.2818
    ##   524     35       2   0.2662
    ##   533     33       1   0.2583
    ##   550     31       1   0.2503
    ##   558     29       1   0.2419
    ##   567     27       1   0.2334
    ##   574     26       1   0.2247
    ##   583     25       1   0.2160
    ##   613     23       1   0.2068
    ##   624     22       1   0.1976
    ##   641     21       1   0.1884
    ##   643     20       1   0.1791
    ##   654     19       1   0.1698
    ##   655     18       1   0.1605
    ##   687     17       1   0.1512
    ##   689     16       1   0.1420
    ##   705     15       1   0.1328
    ##   707     14       1   0.1236
    ##   728     13       1   0.1145
    ##   731     12       1   0.1055
    ##   735     11       1   0.0967
    ##   765      9       1   0.0876
    ##   791      8       1   0.0787
    ##   814      6       1   0.0684
    ##   883      3       1   0.0530

``` r
jum<-survfit(coxph(Surv(time,status)~age+sex))
```

Problems so far:

-   Factors are not supported
-   No confidence intervals
