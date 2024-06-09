# calcualte the calibration weights
EW.p.M <- function(df,varX,gf,M){
  require(Rsolnp)
  na <- length(unique(df$A))
  fn <- function(x){
    sum(x*log(x))
  }
  constraints <- c(M,1)
  p <- NULL
  for (a in 1:na){
    dfa <- df[df$A == a,]
    dfa_x <- dfa[,varX]
    gx_a <- apply(dfa_x, 1, gf)
    eqn <- function(x){
      c <- NULL
      for (i in 1:nrow(gx_a)){
        c <- c(c, sum(x*gx_a[i,]))
      }
      c <- c(c,sum(x))
      return(c)
    }
    nr <- nrow(dfa)
    x0 <- rep(1/nr,nr)
    sol <- solnp(pars = x0, fun = fn, eqfun = eqn, eqB = constraints, 
                 LB = rep(0,nr), UB = rep(1,nr), control = list(trace=0, tol=1e-6))
    p <- c(p, sol$pars)
  }
  
  return(p)
}





# weighted Kaplan Meier RMST estimator
W.KM.Est <- function(df,p,tau){
  EW.mu <- EW.sd <- NULL
  df$A <- as.numeric(df$A) - 1
  mu1 <- my_akm_rmst(df$Y[df$A==1], df$status[df$A==1],p[df$A==1],tau)
  mu0 <- my_akm_rmst(df$Y[df$A==0], df$status[df$A==0],p[df$A==0],tau)
  delta <- mu0$mu - mu1$mu
  EW.mu <- c(mu0$mu, mu1$mu, delta)
  delta.sd <- sqrt(mu1$V + mu0$V)
  EW.sd <- c(sqrt(mu0$V), sqrt(mu1$V), delta.sd)
  return(list(mu = EW.mu, sd = EW.sd))
}



# weighted G-Formula RMST estimator
W.GF.Est <- function(df, p, varX, gf, tau){
  EW.mu <- EW.sd <- NULL
  df$A <- as.numeric(df$A) - 1
  fit <- fit.rmst.reg(df, varX, gf, tau)
  gamma <- coef(fit)
  vgamma <- vcov(fit)
  dfx <- t(apply(df[,varX], 1, gf))
  mdf0 <- cbind(1,0,dfx, 0*dfx)
  mdf1 <- cbind(1,1,dfx, 1*dfx)
  mu0 <- sum(p*(as.matrix(mdf0)%*%matrix(gamma,ncol = 1)))/sum(p)
  mu1 <- sum(p*(as.matrix(mdf1)%*%matrix(gamma,ncol = 1)))/sum(p)
  delta <- mu0 - mu1
  J0 <- t(as.matrix(mdf0)) %*% p/sum(p)
  J1 <- t(as.matrix(mdf1)) %*% p/sum(p)
  delta.J <- t(as.matrix(mdf0)-as.matrix(mdf1)) %*% p/sum(p)
  sd0 <- as.numeric(sqrt(t(J0) %*% vgamma %*% J0))
  sd1 <- as.numeric(sqrt(t(J1) %*% vgamma %*% J1))
  delta.sd <- as.numeric(sqrt(t(delta.J) %*% vgamma %*% delta.J))  
  EW.mu <- c(mu0,mu1,delta)
  EW.sd <- c(sd0,sd1,delta.sd)
  
  return(list(mu = EW.mu, sd = EW.sd))
}


# weighted Hajek RMST estimator
W.HJ.Est <- function(df, p, tau){
  require(geex)
  m_fun <- function(data){
    p <- data$p
    A <- data$A
    Y <- data$Y
    w <- data$w
    function(theta){
      c(p*A*w*(Y-theta[1]),
        p*(1-A)*w*(Y-theta[2]))
    }
  }
  df$A <- as.numeric(df$A) - 1
  A <- df$A
  y <- pmin(df$Y, tau)
  d <- df$status
  d[y==tau]=1
    
  d1=d[A==1]; d0=d[A==0]
  y1=y[A==1]; y0=y[A==0]
    
  fit1=my.func_surv(y1, 1-d1)
  fit0=my.func_surv(y0, 1-d0)
    
  w1=d1/rep(pmax(fit1$surv,0.001), table(y1))
  w0=d0/rep(pmax(fit0$surv,0.001), table(y0))
    
    
  dat_mod <- data.frame(p = c(p[A==1], p[A==0]),
                        w = c(w1, w0),
                        Y = c(y1, y0),
                        A = c(A[A==1], A[A==0]))
    
    
  res <- geex::m_estimate(
         estFUN = m_fun,
         data = dat_mod,
         root_control = setup_root_control(start = c(0.5,0.5))
  )
    
  coef <- coef(res)
  mu1 <- coef[1]
  mu0 <- coef[2]
  delta <- mu0 - mu1
    
  vcov <- vcov(res)
  sd1 <- sqrt(vcov[1,1])
  sd0 <- sqrt(vcov[2,2])
  c1 = matrix(c(-1,1), nrow = 1, ncol = 2)
  delta.sd <- as.numeric(sqrt(c1%*%vcov%*%t(c1)))
    
  EW.mu <- c(mu0,mu1,delta)
  EW.sd <- c(sd0,sd1,delta.sd)
  
  return(list(mu = EW.mu, sd = EW.sd))
}



# weighted Augmented RMST estimator
W.AG.Est <- function(df, p, varX, gf, tau){
  require(geex)
  m_fun <- function(data){
    p <- data$p
    A <- data$A
    Y <- data$Y
    w <- data$w
    mu1 <- data$mu1
    mu0 <- data$mu0
    function(theta){
      c(p*A*w*(Y-mu1-theta[1]),
        p*(1-A)*w*(Y-mu0-theta[2]),
        p*(mu1-theta[3]),
        p*(mu0-theta[4]))
    }
  }
  df$A <- as.numeric(df$A) - 1
  A <- df$A
  y <- pmin(df$Y, tau)
  d <- df$status
  d[y==tau]=1
  
  d1=d[A==1]; d0=d[A==0]
  y1=y[A==1]; y0=y[A==0]
  
  fit1=my.func_surv(y1, 1-d1)
  fit0=my.func_surv(y0, 1-d0)
  
  w1=d1/rep(pmax(fit1$surv,0.001), table(y1))
  w0=d0/rep(pmax(fit0$surv,0.001), table(y0))
  
  
  fit <- fit.rmst.reg(df, varX, gf, tau)
  gamma <- coef(fit)
  vgamma <- vcov(fit)
  dfx <- t(apply(df[,varX], 1, gf))
  mdf0 <- cbind(1,0,dfx, 0*dfx)
  mdf1 <- cbind(1,1,dfx, 1*dfx)
  #mdf0 <- cbind(1,0,df[,paste0("X",1:nX)], 0*df[,paste0("X",1:nX)])
  #mdf1 <- cbind(1,1,df[,paste0("X",1:nX)], 1*df[,paste0("X",1:nX)])
  m0 <- as.matrix(mdf0)%*%matrix(gamma,ncol = 1); m0 <- as.vector(m0)
  m1 <- as.matrix(mdf1)%*%matrix(gamma,ncol = 1); m1 <- as.vector(m1)
  

  dat_mod <- data.frame(p = c(p[A==1], p[A==0]),
                        w = c(w1, w0),
                        Y = c(y1, y0),
                        A = c(A[A==1], A[A==0]),
                        mu1 = c(m1[A==1], m1[A==0]),
                        mu0 = c(m0[A==1], m0[A==0]))
  
  dat_mod <- data.frame(p = c(p[A==1], p[A==0]),
                        w = c(w1, w0),
                        Y = c(y1, y0),
                        A = c(A[A==1], A[A==0]),
                        mu1 = 8,
                        mu0 = 5)
  
  
  res <- geex::m_estimate(
    estFUN = m_fun,
    data = dat_mod,
    root_control = setup_root_control(start = c(1,1,0,0))
  )
  
  coef <- coef(res)
  mu1 <- coef[1] + coef[3]
  mu0 <- coef[2] + coef[4]
  delta <- mu0 - mu1
  
  vcov <- vcov(res)
  c1 <- matrix(c(1,1), nrow = 1, ncol = 2)
  c2 <- matrix(c(-1,1,-1,1), nrow = 1, ncol = 4)
  sd1 <- as.numeric(sqrt(c1%*%vcov[c(1,3),c(1,3)]%*%t(c1)))
  sd0 <- as.numeric(sqrt(c1%*%vcov[c(2,4),c(2,4)]%*%t(c1)))
  delta.sd <- as.numeric(sqrt(c2%*%vcov%*%t(c2)))
  
  EW.mu <- c(mu0,mu1,delta)
  EW.sd <- c(sd0,sd1,delta.sd)
  
  return(list(mu = EW.mu, sd = EW.sd))
}




#####################################################################################
# other source functions
my.func_surv <- function(y, d){
  #--input--
  #y=time
  #d=status
  
  #--
  id=order(y)
  y=y[id]
  d=d[id]
  
  #--
  t_idx = unique(c(0,y))
  ny = length(y)
  
  #--
  Y = N = C = S = H = D = E = rep(0,length(t_idx))
  
  #i=1
  Y[1] = ny
  N[1] = 0
  C[1] = 0
  S[1] = 1
  H[1] = 0
  D[1] = 0
  E[1] = 0
  
  #i>=2
  for(i in 2:length(t_idx)){
    Y[i] = Y[i-1] - N[i-1] - C[i-1]
    N[i] = ifelse(sum(y==t_idx[i] & d==1)>0, sum(y==t_idx[i] & d==1), 0)
    C[i] = ifelse(sum(y==t_idx[i] & d==0)>0, sum(y==t_idx[i] & d==0), 0)
    
    if(Y[i]<0){Y[i] = 0}
    
    S[i] = ifelse(Y[i]==0, S[i-1], S[i-1]*(1-(N[i]/Y[i])))
    H[i] = ifelse(Y[i]*(Y[i]-N[i])==0, 0, N[i]/(Y[i]*(Y[i]-N[i])))
    
    if(S[i]<0){S[i] = 0}
    
    D[i] = sum(H[2:i])
    E[i] = sqrt((S[i]**2)*D[i])
    
    if(is.na(S[i])){S[i] = 0}
    if(is.na(E[i])){E[i] = 0}
  }
  
  #--output--
  out           = as.data.frame(cbind(t_idx, Y, N, C, S, E))
  colnames(out) = c("t_idx", "n_risk", "n_event", "n_censor", "surv", "se")
  
  #--to match the output of survfit--
  out2 = out[t_idx!=0,]
  
  #--
  Z2 = list()
  Z2$out      = out2
  Z2$t_idx    = out2[,"t_idx"]
  Z2$n_risk   = out2[,"n_risk"]
  Z2$n_event  = out2[,"n_event"]
  Z2$n_censor = out2[,"n_censor"]
  Z2$surv     = out2[,"surv"]
  Z2$se       = out2[,"se"]
  
  return(Z2)
}



my.rmst2reg=function(y, delta, arm, x, tau, w=rep(1,length(y))){
  
  n=length(y)
  x=as.matrix(cbind(1, x))
  p=length(x[1,])
  
  y0=pmin(y, tau)
  d0=delta
  d0[y0==tau]=1
  
  d10=d0[arm==1]
  d00=d0[arm==0]
  y10=y0[arm==1]
  y00=y0[arm==0]
  x1=x[arm==1,]
  x0=x[arm==0,]
  n1=length(d10)
  n0=length(d00)
  
  id1=order(y10)
  y10=y10[id1]
  d10=d10[id1]
  x1=x1[id1,]
  
  id0=order(y00)
  y00=y00[id0]
  d00=d00[id0]
  x0=x0[id0,]
  
  fitc1=my.func_surv(y10, 1-d10)
  fitc0=my.func_surv(y00, 1-d00)
  
  weights1=d10/rep(pmax(fitc1$surv,0.001), table(y10))
  weights0=d00/rep(pmax(fitc0$surv,0.001), table(y00))
  
  w1=w[arm==1]
  w0=w[arm==0]
  w1=w1[id1]
  w0=w0[id0]
  weights=c(weights1, weights0)*c(w1,w0)
  
  
  fitt=lm(c(y10,y00)~ rbind(x1, x0)-1, weights=weights)
  
  return(fitt)
}


fit.rmst.reg <- function(df, varX, gf, tau){
  dfx <- t(apply(df[,varX], 1, gf))
  cov <- as.data.frame(cbind(df$A, dfx, df$A*dfx))
  rmst_fit <- my.rmst2reg(y = df$Y,
                          delta = df$status,
                          x = cov,
                          arm = df$A,
                          tau = tau)
  return(rmst_fit)
}




my_akm_rmst <- function(time, status, weight=NULL, tau=NULL){
  
  data <- data.frame(time, status, weight)
  
  #--- AKM ---
  # Based on 'adjusted.KM' function from {IPWsurvival} package
  # Author: F. Le Borgne and Y. Foucher
  tj <- c(0,sort(unique(data$time[data$status==1])))
  dj <- sapply(tj, function(x){sum(data$weight[data$time==x & data$status==1])})
  yj <- sapply(tj, function(x){sum(data$weight[data$time>=x])})
  st <- cumprod(1-(dj/yj))
  m <- sapply(tj, function(x){sum((data$weight[data$time>=x])^2)})
  mj <- ((yj^2)/m)
  #ft <- data.frame(time=tj, n_risk=yj, n_event=dj, survival=st, variable=i, m=mj)
  ft <- data.frame(tj, yj, dj, st, mj)
  
  #--- RMST ---
  # Based on 'rmst1 function' from {survRM2} package
  # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
  rtime <- ft$tj<=tau
  tj_r <- sort(c(ft$tj[rtime],tau))
  st_r <- ft$st[rtime]
  yj_r <- ft$yj[rtime]
  dj_r <- ft$dj[rtime]
  time_diff <- diff(c(0, tj_r))
  areas <- time_diff * c(1, st_r)
  rmst <- sum(areas)
  
  #--- Variance ---
  mj_r <- ft$mj[rtime]
  var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
  #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
  var_r <- c(var_r,0)
  rmst_var <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
  
  return(data.frame(mu = rmst, V = rmst_var))     
}






