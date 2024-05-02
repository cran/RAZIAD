#' Regression Analysis for Zero-altered or Zero-inflated Data
#'
#' @param x a design matrix containing an intercept column (all ones) along with other available covariates available and selected for the response variable
#' @param y a zero-inflated or zero-altered(Hurdle) count response variable, represented as an integer vector
#' @param b0 the initial parameters for the model, calculated as the product of the number of parameters in the specified models and the number of covariates. For simplicity one may put the MLE and intercept of the parameters and set the rest of covariates to zero or change them.
#' @param dist can be specified as follows: "ZIP" for "zero-inflated Poisson",
#' "ZINB" for "zero-inflated negative binomial",
#' "ZINB-r" for "zero-inflated negative binomial with fixed r",
#' "ZIBNB" for "zero-inflated beta negative binomial",
#' "ZIBB" for "zero-inflated beta binomial",
#' "ZIBB-n" for "zero-inflated beta binomial with fixed n",
#' "ZIBB-ab" for "zero-inflated beta binomial with fixed alpha and beta",
#' "PH" for "zero-altered(hurdle) Poisson",
#' "NBH" for "zero-altered(Hurdle) negative binomial",
#' "NBH-r" for "zero-altered (Hurdle) negative binomial with fixed r",
#' "BNBH" for " zero-altered (Hurdle) beta negative binomial",
#' "BBH" for "zero-altered (Hurdle) beta binomial",
#' "BBH-n" for "zero-altered(Hurdle) beta binomial with fixed n",
#'  and "BBH-ab" for "zero-altered(Hurdle) beta binomial with fixed alpha and beta".
#' @param link can be set to one of four different options: "logit" for the logistic link function,
#' "probit" for the probit link function,
#' "loglog" for the log-log link function,
#' and "cloglog" for the complementary log-log link function
#' @return  A list containing AIC, BIC, the corresponding value of log likelihood, and the maximum likelihood estimate (MLE) of the unknown parameters in the model.
#'  If dist = ZIP, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log.likelihhod: The value of log likelihood with maximum likelihood estimate plugged-in.
#' \item Estimated.Parameters: The maximum likelihood estimate of \eqn{\Gamma} and \eqn{\beta} for intercept and covariates included in the design matrix.
#' }
#' If dist = ZINB, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\Gamma}), (\eqn{\beta_1}), and (\eqn{\beta_2}) for intercept and covariates included in the design matrix.
#' }
#' If dist = ZINB-r, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\beta_1}) for intercept only and (\eqn{\Gamma}), and (\eqn{\beta_2}) for intercept and covariates included in the design matrix.
#' }
#' If dist = ZIBNB, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\Gamma}), (\eqn{\beta_1}), (\eqn{\beta_2}), and (\eqn{\beta_3}) for intercept and covariates included in the design matrix.
#' }
#' If dist = ZIBB, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\Gamma}), (\eqn{\beta_1}), (\eqn{\beta_2}), and (\eqn{\beta_3}) for intercept and covariates included in the design matrix.
#' }
#' If dist = ZIBB-n, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\beta_1}) for intercept only and (\eqn{\Gamma}),  (\eqn{\beta_2}), and (\eqn{\beta_3}) for intercept and covariates included in the design matrix.
#' }
#' If dist = ZIBB-ab, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\beta_2}) and (\eqn{\beta_3}) for intercept only and (\eqn{\Gamma}) and (\eqn{\beta_1}) for intercept and covariates included in the design matrix.
#' }
#' If dist = PH, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item  Estimated.Parameters: The maximum likelihood estimate of \eqn{\Gamma} and \eqn{\beta} for intercept and covariates included in the design matrix.
#' }
#' If dist = NBH, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\Gamma}), (\eqn{\beta_1}), and (\eqn{\beta_2}) for intercept and covariates included in the design matrix.
#' }
#' If dist = NBH-r, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\beta_1}) for intercpet only anf (\eqn{\Gamma}), and (\eqn{\beta_2}) for intercept and covariates included in the design matrix.
#' }
#' If dist = BNBH, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\Gamma}), (\eqn{\beta_1}), (\eqn{\beta_2}), and (\eqn{\beta_3}) for intercept and covariates included in the design matrix.
#' }
#' If dist = BBH, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\Gamma}), (\eqn{\beta_1}), (\eqn{\beta_2}), and (\eqn{\beta_3}) for intercept and covariates included in the design matrix.
#' }
#' If dist = BBH-n, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\beta_1}) for intercpet only and (\eqn{\Gamma}), (\eqn{\beta_2}), and (\eqn{\beta_3}) for intercept and covariates included in the design matrix.
#' }
#' If dist = BBH-ab, the following values are returned:
#' \itemize{
#' \item AIC: Akaike Information Criterion, a measure of the model's goodness of fit adjusted for the number of parameters.
#' \item BIC: Bayesian Information Criterion, a criterion for model selection among a finite set of models.
#' \item log-likelihood: The value of log-likelihood corresponding to the maximum likelihood estimate.
#' \item Estimated Parameters: Maximum likelihood estimates of (\eqn{\beta_2}) and (\eqn{\beta_3}) for intercept only and (\eqn{\Gamma}) and (\eqn{\beta_1}), for intercept and covariates included in the design matrix.
#' }
#' @export
#' @examples
#' intercept<- rep(1,4406)
#'dt = DebTrivedi[, c(6:8, 13, 15, 18)]
#'dt = cbind(intercept, dt)
#'dt$gender.male <- ifelse(dt$gender == 'male', 1, 0)
#'dt$gender.female <- ifelse(dt$gender == 'female', 1, 0)
#'dt$health.poor <- ifelse(dt$health == 'poor', 1, 0)
#'dt$health.average <- ifelse(dt$health == 'average', 1, 0)
#'dt$health.excellent <- ifelse(dt$health == 'excellent', 1, 0)
#'dt$privins.yes <- ifelse(dt$privins == 'yes', 1, 0)
#'dt$privins.no <- ifelse(dt$privins == 'no', 1, 0)
#'y = DebTrivedi[,1]
#'x = data.matrix(dt[, c(1, 2, 4, 6, 8, 10, 12, 13)])
#'np = dim(x)[2]
#'b0 = c(rep(0.3, np), rep(0.1, np))
#'reg.model(x, y, b0=b0, dist="PH", link="logit")
reg.model <- function(x, y, b0=NULL, dist = "ZIP", link)
  {
  if (link == "loglog")  {
    ginv = function(phi) {exp(-exp(phi));};                       # g-inverse with "log-log" link
    ginvd = function(phi){-exp(phi-exp(phi));};                   # g-inverse derivative
  }
  if (link == "cloglog") {
    ginv = function(phi) {1-exp(-exp(phi));};                     # g-inverse with "c-log-log"
    ginvd = function(phi){exp(phi-exp(phi));};                    # g-inverse derivative
  }
  if (link == "logit")   {
    ginv = function(phi){exp(phi)/(1 + exp(phi));};               # g-inverse with "logit" link
    ginvd = function(phi){exp(phi)/(1+exp(phi))^2;};              # g-inverse derivative
  }
  if (link == "probit")  {
    ginv = function(phi){stats::pnorm(phi, mean = 0, sd = 1);};   # g-inverse with "probit" link
    ginvd = function(phi){stats::dnorm(phi, mean = 0, sd = 1);};  # g-inverse derivative
  }
  if (dist == "PH")      {
    hinv = function(theta){exp(theta);};                                     # h-inverse
    hinvd = function(theta) {exp(theta);};                                   # h-inverse derivative
    p0 = function(x){exp(-x);};                                              # P_0(lambda)
    p0d = function(x){-exp(-x);};                                            # P_0(lambda) derivative
    ptheta = function(x){stats::dpois(x[1], x[2], log=F);};
    pthetad = function(x){stats::dpois(x[1], x[2], log=F)*((x[1]/x[2])-1);};
    np = dim(x)[2];             # number of predictors
    N = dim(x)[1];
    if(is.null(b0)) {b0 = rep(0.1, 2*np)}
    {
      neg.log.lik <- function(theta) {
        g <- theta[1:np]        # gamma
        b <- theta[np + (1:np)] # beta
        xg = x%*%g;             # X^T*gamma
        xb = x%*%b;             # X^T*beta
        phii = ginv(xg);        # phi_i
        thetai = hinv(xb);      # theta_i
        thetay = cbind(y, thetai);
        thetay = thetay[y>0,];
        ans = -sum(log(phii)[y == 0]) - sum(log(1 - phii)[y>0]) -
               sum(log(apply(thetay,1,ptheta) + 1e-15)) + sum (log(1 - p0(thetai))[y>0])
        return(ans)
      }
      gp <- function(theta){
        g = theta[1:np]         # gamma
        b = theta[np + (1:np)]  # beta
        xg = x%*%g;             # X^T*gamma
        xb = x%*%b;             # X^T*beta
        phii = ginv(xg);        # phi_i
        phiid = ginvd(xg)       # phi_i derivative
        thetai = hinv(xb);      # theta_i
        thetaid = hinvd(xb);    # theta_i derivative
        thetay = cbind(y, thetai);
        dl = - apply((phiid/(phii))[y==0] * x[y==0,], 2, sum) + apply((phiid/(1-phii))[y>0] * x[y>0,], 2, sum)
        df = -apply((apply(thetay, 1, pthetad) * thetaid/(apply(thetay, 1, ptheta)+1e-15))[y>0] * x[y>0,], 2, sum) -
          apply (((p0d(thetai) * thetaid)/(1 - p0(thetai)))[y>0] * x[y>0,], 2, sum)
        return(c(dl , df))
      }
      estimate = stats::optim(par = b0,fn = neg.log.lik, gr = gp, method = "BFGS")
      estlike = estimate$value
      mle = matrix(c(estimate$par), nrow = 1)
      Akaike.information.criterion = -2* -(estlike) + 2 * 2 * np
      Bayesian.information.criterion = -2 * -(estlike) + log (N) * 2 * np
      vec1 <- 1:np
      gamma <- mle[vec1]
      vec2 <- np+1:np
      beta <- mle[vec2]
      Estimated = matrix (c(gamma, beta), ncol = np, nrow = 2, byrow = TRUE)
      rownames(Estimated) = c("Estimated Gamma", "Estimated Beta")
      return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike , Estimated.Parameters = Estimated ))
    }
  }
  if (dist == "ZIP")     {
    hinv = function(theta) {exp(theta);};                                 # h-inverse
    hinvd = function(theta) {exp(theta);};                                # h-inverse derivative
    p0 = function(x){exp(-x);};                                 # P_0(lambda)
    p0d = function(x){-exp(-x);};                               # P_0(lambda) derivative
    ptheta = function(x){stats::dpois(x[1], x[2]);};                      # c(theta, y)
    pthetad = function(x){stats::dpois(x[1], x[2])*((x[1]/x[2])-1);};     # c(theta,y) derivative
    np = dim(x)[2];              # number of predictors
    N = dim(x)[1];
    if(is.null(b0)) {b0 = rep(0.1,2*np)}
    neg.log.lik <- function(theta) {
      g <- theta[1:np]           # gamma
      b <- theta[np+(1:np)]      # beta
      xg = x%*%g;                # X^T*gamma
      xb = x%*%b;                # X^T*beta
      phii = ginv(xg);           # phi_i
      thetai = hinv(xb);         # theta_i
      thetay = cbind(y, thetai);
      thetay = thetay[y>0,];
      ans = -sum(log(1 - phii)[y>0]) - sum(log(apply(thetay, 1, ptheta) + 1e-15)) -
             sum(log(phii + (1 - phii) * p0(thetai))[y==0])
      return(ans)
    }
    gp <- function(theta) {
      g <- theta[1:np];          # gamma
      b <- theta[np+(1:np)];     # beta
      xg = x%*%g;                # X^T*gamma
      xb = x%*%b;                # X^T*beta
      phii = ginv(xg);           # phi_i
      phiid = ginvd(xg)          # phi_i derivative
      thetai = hinv(xb);         # theta_i
      thetaid=hinvd(xb);         # theta_i derivative
      thetay = cbind(y, thetai);
      dl = apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum) -
        apply(((1 - p0(thetai)) * phiid/(phii + (1 - phii) * p0(thetai)))[y==0] * x[y==0,], 2, sum)
      df = -apply((apply(thetay, 1, pthetad) * thetaid/(apply(thetay, 1, ptheta) + 1e-15))[y>0] * x[y>0,], 2, sum) -
        apply((((1 - phii) * p0d(thetai) * thetaid)/(phii + (1 - phii) * p0(thetai)))[y==0] * x[y==0,], 2, sum)
      return(c(dl, df))
    }
    estimate = stats::optim(par = b0,fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    mle = matrix(c(estimate$par), nrow = 1)
    Akaike.information.criterion = -2* -(estlike) + 2 * 2 * np
    Bayesian.information.criterion = -2 * -(estlike) + log (N) * 2 * np
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta <- mle[vec2]
    Estimated = matrix (c(gamma, beta), ncol = np, nrow = 2, byrow = TRUE)
    rownames(Estimated) = c("Estimated Gamma", "Estimated Beta")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike , Estimated.Parameters = Estimated ))
  }
  if (dist == "NBH")     {
    h1inv = function(theta) {exp(theta);};                                      # h1-inverse
    h1invd = function(theta) {exp(theta);};                                     # h1-inverse derivative
    h2inv = function(theta) {exp(theta)/(1+exp(theta));};                       # h2-inverse
    h2invd = function(theta){exp(theta)/(1+exp(theta))^2;};                     # h2-inverse derivative
    p01 = function(x) {x[2]^x[1];}; #x[2]=p                                      # P_0(r,p)
    p0d.r = function(x) {log(x[2]) * (x[2]^x[1])}                               # P_0(r) derivative
    p0d.p = function(x) {x[1] * ((x[2])^(x[1]-1))}                              # P_0(p) derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size = y[2], prob= y[3]);};       # f(theta, y)
    #ptheta=function(y){dnbinom(y[1], size= y[2], prob=y[3]+1e-20);};
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3]) * (digamma(x[1] + x[2]) + log(x[3]) - digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3]) * (x[2]/x[3] - (x[1]/(1-x[3])));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    if(is.null(b0)) {b0=rep(0.2,3*np);}
    neg.log.lik <- function(theta) {
      g <- theta[1:np];          # gamma
      b1 <- theta[np+(1:np)];    # beta1
      b2 <- theta[2*np+(1:np)];  # beta2
      xg = x%*%g;                # X^T*gamma
      xb1 = x%*%b1;              # X^T*beta1
      xb2 = x%*%b2;              # X^T*beta2
      phii = ginv(xg);           # phi_i
      theta1i = h1inv(xb1)       # theta_1 i
      theta2i = h2inv(xb2)       # theta_2 i
      thetai0 = cbind(theta1i, theta2i)
      thetay = cbind(y, theta1i, theta2i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #loglik = -Inf;
      #if((min(thetai0[,2])>0)&(max(thetai0[,2])<1)&(min(thetai0[,1])>0)&(min(phii)>0)&(max(phii)<1)) {
      ans = -sum(log(phii)[y==0]) - sum(log(1-phii)[y>0]) -
        sum(log(apply(thetay1, 1, ptheta1))) + sum(log(1 - apply(thetai01, 1, p01)));
      return(ans)
      #}
    }
    gp <- function(theta) {
      g <- theta[1:np];           # gamma
      b1 <- theta[np+(1:np)];     # beta1
      b2 <- theta[2*np+(1:np)];   # beta2
      xg = x%*%g;                 # X^T*gamma
      xb1 = x%*%b1;               # X^T*beta1
      xb2 = x%*%b2;               # X^T*beta2
      phii = ginv(xg);            # g-inverse
      phiid = ginvd(xg)           # g-inverse derivative
      theta1i = h1inv(xb1)        # h1-inverse
      theta1id = h1invd(xb1)      # h1-inverse derivative
      theta2i = h2inv(xb2)        # h2-inverse
      theta2id = h2invd(xb2)      # h2-inverse derivative
      thetai0 = cbind(theta1i, theta2i)
      thetay = cbind(y, theta1i, theta2i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,]
      thetai01 = thetai0[y>0,];
      #wrt phi
      dl = -apply((phiid/(phii))[y==0] * x[y==0,], 2, sum) + apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      #wrt r
      c1a = apply(thetay1, 1, pthetad.r) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta1)
      c1c = (c1a/c1b) * x[y>0,]
      temp3.r = apply(c1c, 2, sum)
      d1a = apply(thetai01, 1, p0d.r) * (theta1id[y>0])
      d1b = 1 - (apply(thetai01, 1 ,p01))
      d1c = (d1a/d1b) * x[y>0,]
      temp4.r = apply(d1c, 2, sum)
      dr = -temp3.r - temp4.r
      #wrt p
      e1 = apply(thetay1, 1, pthetad.p) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta1)
      e3 = (e1/e2) * x[y>0,]
      temp3.p = apply(e3, 2, sum)
      f1 = apply(thetai01, 1, p0d.p) * theta2id[y>0]
      f2 = 1 - (apply(thetai01, 1, p01))
      f3 = (f1/f2) * x[y>0,]
      temp4.p = apply(f3, 2, sum)
      dp = -temp3.p - temp4.p
      return(c(dl, dr, dp))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    mle = matrix(c(estimate$par), nrow = 1)
    Akaike.information.criterion = -2* -(estlike) + 2 * 3 * np
    Bayesian.information.criterion = -2 * -(estlike) + log (N) * 3 * np
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta1 <- mle[vec2]
    vec3 <- 2*np+(1:np)
    beta2 <- mle[vec3]
    Estimated = matrix (c(gamma, beta1, beta2), ncol = np, nrow = 3, byrow = TRUE)
    rownames(Estimated) = c("Estimated Gamma", "Estimated Beta1", "Estimated Beta2")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike , Estimated.Parameters = Estimated ))
  }
  if (dist == "NBH-r")   {
    h1inv = function(theta) {exp(theta);};                                      # h1-inverse
    h1invd = function(theta) {exp(theta);};                                     # h1-inverse derivative
    h2inv = function(theta) {exp(theta)/(1+exp(theta));};                       # h2-inverse
    h2invd = function(theta){exp(theta)/(1+exp(theta))^2;};                     # h2-inverse derivative
    p01 = function(x) {x[2]^x[1];}; #x[2]=p                                      # P_0(r,p)
    p0d.r = function(x) {log(x[2]) * (x[2]^x[1])}                               # P_0(r) derivative
    p0d.p = function(x) {x[1] * ((x[2])^(x[1]-1))}                              # P_0(p) derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size = y[2], prob= y[3]);};      # f(theta, y)
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3]) * (digamma(x[1] + x[2]) + log(x[3]) - digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3]) * (x[2]/x[3] - (x[1]/(1-x[3])));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    xa = as.matrix(rep(1, N))
    if(is.null(b0)) {b0=rep(0.2,3*np);}
    neg.log.lik <- function(theta) {
      g <- theta[1:np];          # gamma
      b1 <- theta[np+1];         # beta1
      b2 <- theta[np+1+(1:np)];  # beta2
      xg = x%*%g;                # X^T*gamma
      xb1 = xa%*%b1;             # X^T*beta1
      xb2 = x%*%b2;              # X^T*beta2
      phii = ginv(xg);           # phi_i
      theta1i = h1inv(xb1)       # theta_1 i
      theta2i = h2inv(xb2)       # theta_2 i
      thetai0 = cbind(theta1i, theta2i)
      thetay = cbind(y, theta1i, theta2i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #loglik = -Inf;
      #if((min(thetai0[,2])>0)&(max(thetai0[,2])<1)&(min(thetai0[,1])>0)&(min(phii)>0)&(max(phii)<1)) {
      ans = -sum(log(phii)[y==0]) - sum(log(1-phii)[y>0]) -
        sum(log(apply(thetay1, 1, ptheta1))) + sum(log(1 - apply(thetai01, 1, p01)));
      return(ans)
      #}
    }
    gp <- function(theta) {
      g <- theta[1:np];           # gamma
      b1 <- theta[np+1];          # beta1
      b2 <- theta[np+1+(1:np)];   # beta2
      xg = x%*%g;                 # X^T*gamma
      xb1 = xa%*%b1;              # X^T*beta1
      xb2 = x%*%b2;               # X^T*beta2
      phii = ginv(xg);            # g-inverse
      phiid = ginvd(xg)           # g-inverse derivative
      theta1i = h1inv(xb1)        # h1-inverse
      theta1id = h1invd(xb1)      # h1-inverse derivative
      theta2i = h2inv(xb2)        # h2-inverse
      theta2id = h2invd(xb2)      # h2-inverse derivative
      thetai0 = cbind(theta1i, theta2i)
      thetay = cbind(y, theta1i, theta2i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,]
      thetai01 = thetai0[y>0,];
      #wrt phi
      dl = -apply((phiid/(phii))[y==0] * x[y==0,], 2, sum) + apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      #wrt r
      c1a = apply(thetay1, 1, pthetad.r) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta1)
      c1c = (c1a/c1b) * xa[y>0,]
      temp3.r = sum(c1c)
      #temp3.r = apply(c1c, 2, sum)
      d1a = apply(thetai01, 1, p0d.r) * (theta1id[y>0])
      d1b = 1 - (apply(thetai01, 1 ,p01))
      d1c = (d1a/d1b) * xa[y>0,]
      temp4.r = sum(d1c)
      #temp4.r = apply(d1c, 2, sum)
      dr = -temp3.r - temp4.r
      #wrt p
      e1 = apply(thetay1, 1, pthetad.p) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta1)
      e3 = (e1/e2) * x[y>0,]
      temp3.p = apply(e3, 2, sum)
      f1 = apply(thetai01, 1, p0d.p) * theta2id[y>0]
      f2 = 1 - (apply(thetai01, 1, p01))
      f3 = (f1/f2) * x[y>0,]
      temp4.p = apply(f3, 2, sum)
      dp = -temp3.p - temp4.p
      return(c(dl, dr, dp))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    mle = matrix(c(estimate$par), nrow = 1)
    Akaike.information.criterion = -2* -(estlike) + 2 * 3 * np
    Bayesian.information.criterion = -2 * -(estlike) + log (N) * 3 * np
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+1
    beta1 <- mle[vec2]
    vec3 <- np+1+(1:np)
    beta2 <- mle[vec3]
    #Estimated = matrix (c(gamma, beta1, beta2), ncol = np, nrow = 3, byrow = TRUE)
    #rownames(Estimated) = c("Estimated Gamma", "Estimated Beta1", "Estimated Beta2")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion,
                log.likelihhod = estlike , Estimated.Gamma = gamma, Estimated.Beta1 = beta1, Estimated.Beta2 = beta2 ))
  }
  if (dist == "ZINB")    {
    h1inv = function(theta) {exp(theta);};                                      # h1-inverse
    h1invd = function(theta) {exp(theta);};                                     # h1-inverse derivative
    h2inv = function(theta) {exp(theta)/(1+exp(theta));};                       # h2-inverse
    h2invd = function(theta){exp(theta)/(1+exp(theta))^2;};                     # h2-inverse derivative
    p01 = function(x) {x[2]^x[1];}; #x[2]=p                                      # P_0(r,p)
    p0d.r = function(x) {log(x[2]) * (x[2]^x[1])}                               # P_0(r) derivative
    p0d.p = function(x) {x[1] * ((x[2])^(x[1]-1))}                              # P_0(p) derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size = y[2], prob= y[3]);};       # f(theta, y)
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3]) * (digamma(x[1] + x[2]) + log(x[3]) - digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3]) * (x[2]/x[3] - (x[1]/(1-x[3])));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    neg.log.lik <- function(theta) {
      g <- theta[1:np];          # gamma
      b1 <- theta[np+(1:np)];    # beta1
      b2 <- theta[2*np+(1:np)];  # beta2
      xg = x%*%g;                # X^T*gamma
      xb1 = x%*%b1;              # X^T*beta1
      xb2 = x%*%b2;              # X^T*beta2
      phii = ginv(xg);           # phi_i
      theta1i = h1inv(xb1)       # theta_1 i
      theta2i = h2inv(xb2)       # theta_2 i
      thetai0 = cbind(theta1i, theta2i)
      thetay = cbind(y, theta1i, theta2i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      ans = -sum(log(1 - phii)[y>0]) - sum(log(apply(thetay1, 1, ptheta1))) -
             sum(log(phii[y==0] + (1 - phii[y==0]) * apply(thetai00, 1, p01)));
      return(ans)
    }
    gp <- function(theta) {
      g <- theta[1:np];           # gamma
      b1 <- theta[np+(1:np)];     # beta1
      b2 <- theta[2*np+(1:np)];   # beta2
      xg = x%*%g;                 # X^T*gamma
      xb1 = x%*%b1;               # X^T*beta1
      xb2 = x%*%b2;               # X^T*beta2
      phii = ginv(xg);            # g-inverse
      phiid = ginvd(xg)           # g-inverse derivative
      theta1i = h1inv(xb1)        # h1-inverse
      theta1id = h1invd(xb1)      # h1-inverse derivative
      theta2i = h2inv(xb2)        # h2-inverse
      theta2id = h2invd(xb2)      # h2-inverse derivative
      thetai0 = cbind(theta1i, theta2i)
      thetay = cbind(y, theta1i, theta2i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,]
      thetai01 = thetai0[y>0,];
      #wrt phi
      temp1 = apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      temp2a = (1 - apply(thetai00, 1, p01)) * phiid[y==0]
      temp2b = phii[y==0] + ((1 - phii)[y==0]) * apply(thetai00, 1, p01)
      temp2c = (temp2a/temp2b) * x[y==0,]
      temp2 = apply(temp2c, 2, sum)
      dl = temp1 - temp2
      #wrt r
      c1a = apply(thetay1, 1, pthetad.r) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta1)
      c1c = (c1a/c1b) * x[y>0,]
      temp3.r = apply(c1c, 2, sum)
      d1a = (1 - phii[y==0]) * apply(thetai00, 1, p0d.r) * theta1id[y==0]
      d1b = phii[y==0] + (1 - phii[y==0])*(apply(thetai00, 1, p01))
      d1c = (d1a/d1b) * x[y==0,]
      temp4.r = apply(d1c, 2, sum)
      dr = -temp3.r - temp4.r
      #wrt p
      e1 = apply(thetay1, 1, pthetad.p) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta1)
      e3 = (e1/e2) * x[y>0,]
      temp3.p = apply(e3, 2, sum)
      f1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.p) * theta2id[y==0]
      f2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p01))
      f3 = (f1/f2) * x[y==0,]
      temp4.p = apply(f3, 2, sum)
      dp = -temp3.p - temp4.p
      return(c(dl, dr, dp))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    mle = matrix(c(estimate$par), nrow = 1)
    Akaike.information.criterion = -2* -(estlike) + 2 * 3 * np
    Bayesian.information.criterion = -2 * -(estlike) + log (N) * 3 * np
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta1 <- mle[vec2]
    vec3 <- 2*np+(1:np)
    beta2 <- mle[vec3]
    Estimated = matrix (c(gamma, beta1, beta2), ncol = np, nrow = 3, byrow = TRUE)
    rownames(Estimated) = c("Estimated Gamma", "Estimated Beta1", "Estimated Beta2")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike , Estimated.Parameters = Estimated ))
  }
  if (dist == "ZINB-r")  {
    h1inv = function(theta) {exp(theta);};                                      # h1-inverse
    h1invd = function(theta) {exp(theta);};                                     # h1-inverse derivative
    h2inv = function(theta) {exp(theta)/(1+exp(theta));};                       # h2-inverse
    h2invd = function(theta){exp(theta)/(1+exp(theta))^2;};                     # h2-inverse derivative
    p01 = function(x) {x[2]^x[1];}; #x[2]=p                                      # P_0(r,p)
    p0d.r = function(x) {log(x[2]) * (x[2]^x[1])}                               # P_0(r) derivative
    p0d.p = function(x) {x[1] * ((x[2])^(x[1]-1))}                              # P_0(p) derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size = y[2], prob = y[3]);};     # f(theta, y)
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3]) * (digamma(x[1] + x[2]) + log(x[3]) - digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3]) * (x[2]/x[3] - (x[1]/(1-x[3])));};
    np = dim(x)[2];               # number of predictors
    N = dim(x)[1];
    xa = as.matrix(rep(1, N))
    if(is.null(b0)) {b0=rep(0.2,2*np+1);}
    neg.log.lik <- function(theta) {
      g <- theta[1:np];          # gamma
      b1 <- theta[np+1];         # beta1
      b2 <- theta[np+1+(1:np)];  # beta2
      xg = x%*%g;                # X^T*gamma
      xb1 = xa%*%b1;             # X^T*beta
      xb2 = x%*%b2;              # X^T*beta
      phii = ginv(xg);           # phi_i
      theta1i = h1inv(xb1)       # theta_1 i
      theta2i = h2inv(xb2)       # theta_2 i
      thetai0 = cbind(theta1i, theta2i)
      thetay = cbind(y, theta1i, theta2i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      #    loglik=-Inf;
      #    if((min(thetai0[,2])>0)&(max(thetai0[,2])<1)&(min(thetai0[,1])>0)&(min(phii)>0)&(max(phii)<1)) {
      ans = -sum(log(1 - phii)[y>0]) - sum(log(apply(thetay1, 1, ptheta1))) - sum(log(phii[y==0] + (1 - phii[y==0]) * apply(thetai00, 1, p01)));
      #   };
      return(ans)
    }
    gp <- function(theta) {
      g <- theta[1:np];           # gamma
      b1 <- theta[np+1];          # beta1
      b2 <- theta[np+1+(1:np)];   # beta2
      xg = x%*%g;                 # X^T*gamma
      xb1 = xa%*%b1;              #X^T*beta1
      xb2 = x%*%b2;               #X^T*beta2
      phii = ginv(xg);            # g-inverse
      phiid = ginvd(xg)           # g-inverse derivative
      theta1i = h1inv(xb1)        # h1-inverse
      theta1id = h1invd(xb1)      # h1-inverse derivative
      theta2i = h2inv(xb2)        # h2-inverse
      theta2id = h2invd(xb2)      # h2-inverse derivative
      thetai0 = cbind(theta1i, theta2i)
      thetay = cbind(y, theta1i, theta2i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,]
      thetai01 = thetai0[y>0,];
      #wrt phi
      temp1 = apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      temp2a = (1 - apply(thetai00, 1, p01)) * phiid[y==0]
      temp2b = phii[y==0] + ((1 - phii)[y==0]) * apply(thetai00, 1, p01)
      temp2c = (temp2a/temp2b) * x[y==0,]
      temp2 = apply(temp2c, 2, sum)
      dl = temp1 - temp2
      #wrt r
      c1a = apply(thetay1, 1, pthetad.r) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta1)
      c1c = (c1a/c1b) * xa[y>0,]
      temp3.r = sum(c1c)
      #temp3.r=apply(c1c, 2, sum)
      d1a = (1 - phii[y==0]) * apply(thetai00, 1, p0d.r) * theta1id[y==0]
      d1b = phii[y==0] + (1-phii[y==0]) * (apply(thetai00, 1, p01))
      d1c = (d1a/d1b) * xa[y==0,]
      #temp4.r=apply(d1c, 2, sum)
      temp4.r = sum(d1c)
      dr = -temp3.r - temp4.r
      #wrt p
      e1 = apply(thetay1, 1, pthetad.p) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta1)
      e3 = (e1/e2) * x[y>0,]
      temp3.p = apply(e3, 2, sum)
      f1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.p) * theta2id[y==0]
      f2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p01))
      f3 = (f1/f2) * x[y==0,]
      temp4.p = apply(f3, 2, sum)
      dp = -temp3.p - temp4.p
      return(c(dl, dr, dp))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    mle = matrix(c(estimate$par), nrow = 1)
    Akaike.information.criterion = -2* -(estlike) + 2 * 3 * np
    Bayesian.information.criterion = -2 * -(estlike) + log (N) * 3 * np
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- (np+1)
    beta1 <- mle[vec2]
    vec3 <- np+1+(1:np)
    beta2 <- mle[vec3]
    #Estimated = matrix (c(gamma, beta1, beta2), ncol = np, nrow = 3, byrow = TRUE)
    #rownames(Estimated) = c("Estimated Gamma", "Estimated Beta1", "Estimated Beta2")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion,
                log.likelihhod = estlike , Estimated.Gamma = gamma, Estimated.Beta1 = beta1, Estimated.Beta2 = beta2 ))
  }
  if (dist == "ZIBNB")   {
    h1inv = function(theta) {exp(theta);};    # h1-inverse
    h2inv = function(theta) {exp(theta);};    # h2-inverse
    h3inv = function(theta) {exp(theta);};    # h3-inverse
    h1invd = function(theta) {exp(theta);};   # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};   # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};   # h3-inverse derivative
    ptheta2 = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]);}; #c(y, n, a,b)
    #p0=function(x) {gamma(x[2]+x[4])*gamma(x[3]+x[4])/(gamma(x[2]+x[3]+x[4])*gamma(x[4]));}; #x[1]=n
    p02 = function(x){extraDistr::dbnbinom(0, size = x[2], alpha = x[3], beta = x[4]);}; #c(0, n, a,b)
    pthetad.r = function(x){extraDistr::dbnbinom(x[1], size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[1]) + digamma(x[2]+x[3]) - digamma(x[2]) - digamma(x[1]+x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbnbinom(x[1], size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) + digamma(x[3]+x[4]) - digamma(x[1]+x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbnbinom(x[1], size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[1]+x[3]+x[2]+x[4]) - digamma(x[4]));};
    p0d.r = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) + digamma(x[3]+x[4]) - digamma(x[2]+x[3]+x[4]) + digamma(x[3]));};
    p0d.b = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    if(is.null(b0)) {b0=rep(0.3,4*np);}
    neg.log.lik <- function(theta) {
      g  <- theta[1:np];              # gamma
      b1 <- theta[np+(1:np)];         # beta1
      b2 <- theta[2*np+(1:np)];       # beta2
      b3 <- theta[3*np+(1:np)];       # beta3
      xg = x%*%g;                     # X^T*gamma
      xb1 = x%*%b1;                   # X^T*beta1
      xb2 = x%*%b2;                   # X^T*beta2
      xb3 = x%*%b3;                   # X^T*beta3
      phii = ginv(xg);                # g-inverse
      theta1i = h1inv(xb1)            # theta_1 i
      theta2i = h2inv(xb2)            # theta_2 i
      theta3i = h3inv(xb3)            # theta_3 i
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      ans = -sum(log(1 - phii)[y>0]) - sum(log(apply(thetay1, 1, ptheta2))) + sum(log(phii[y==0] + (1 - phii[y==0]) * apply(thetai00, 1, p02)));
      return(ans)
    }
    gp <- function(theta) {
      g  <- theta[1:np];               # gamma
      b1 <- theta[np+(1:np)];          # beta1
      b2 <- theta[2*np+(1:np)];        # beta2
      b3 <- theta[3*np+(1:np)];        # beta3
      xg = x%*%g;                      # X^T*gamma
      xb1 = x%*%b1;                    # X^T*beta1
      xb2 = x%*%b2;                    # X^T*beta2
      xb3 = x%*%b3;                    # X^T*beta3
      phii = ginv(xg);                 # g-inverse
      theta1i =  ceiling(h1inv(xb1))   # h1-inverse
      theta2i = h2inv(xb2)             # h2-inverse
      theta3i = h3inv(xb3)             # h3-inverse
      phiid = ginvd(xg)                # g-inverse derivative
      theta1id = h1invd(xb1)           # h1-inverse derivative
      theta2id = h2invd(xb2)           # h2-inverse derivative
      theta3id = h3invd(xb3)           # h3-inverse derivative
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00=thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #wrt phi
      temp1 = apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      temp2a = (1 - apply(thetai00, 1, p02)) * phiid[y==0]
      temp2b = phii[y==0] + ((1 - phii)[y==0]) * apply(thetai00, 1, p02)
      temp2c = (temp2a/temp2b) * x[y==0,]
      temp2 = apply(temp2c, 2, sum)
      dl = temp1 - temp2
      # wrt r
      c1a = apply(thetay1, 1, pthetad.r) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta2)
      c1c = (c1a/c1b) * x[y>0,]
      temp3.r = apply(c1c, 2, sum)
      d1a = (1 - phii[y==0]) * apply(thetai00, 1, p0d.r) * theta1id[y==0]
      d1b = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p02))
      d1c = (d1a/d1b) * x[y==0,]
      temp4.r = apply(d1c, 2, sum)
      dr = -temp3.r - temp4.r
      # wrt a
      e1 = apply(thetay1, 1, pthetad.a) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta2)
      e3 = (e1/e2) * x[y>0,]
      temp3.a = apply(e3, 2, sum)
      f1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.a) * theta2id[y==0]
      f2 = phii[y==0] + (1-phii[y==0]) * (apply(thetai00,1,p02))
      f3 = (f1/f2) * x[y==0,]
      temp4.a = apply(f3, 2, sum)
      da = -temp3.a - temp4.a
      # wrt b
      z1 = apply(thetay1, 1, pthetad.b) * (theta3id[y>0])
      z2 = apply(thetay1, 1, ptheta2)
      z3 = (z1/z2) * x[y>0,]
      temp3.b = apply(z3, 2, sum)
      k1 = (1-phii[y==0]) * apply(thetai00, 1, p0d.b) * theta3id[y==0]
      k2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p02))
      k3 = (k1/k2) * x[y==0,]
      temp4.b = apply(k3, 2, sum)
      db = -temp3.b - temp4.b
      return(c(dl, dr, da, db))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method =  "BFGS")
    estlike = estimate$value
    Akaike.information.criterion = -2 * -(estlike) + 2 * 4 * np
    Bayesian.information.criterion = -2 * -(estlike) + (log (N)) * 4 * np
    mle = matrix(c(estimate$par), nrow = 1)
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta1 <- mle[vec2]
    vec3 <- 2*np+(1:np)
    beta2 <- mle[vec3]
    vec4 <- 3*np+(1:np)
    beta3 <- mle[vec4]
    Estimated = matrix (c(gamma, beta1, beta2, beta3), ncol = np, nrow = 4, byrow = TRUE)
    rownames(Estimated) = c("Estimated Gamma", "Estimated Beta1", "Estimated Beta2", "Estimated Beta3")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike , Estimated.Parameters = Estimated ))
  }
  if (dist == "BNBH")    {
    h1inv = function(theta) {exp(theta);};    # h1-inverse
    h2inv = function(theta) {exp(theta);};    # h2-inverse
    h3inv = function(theta) {exp(theta);};    # h3-inverse
    h1invd = function(theta) {exp(theta);};   # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};   # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};   # h3-inverse derivative
    ptheta2 = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]);}; #c(y, n, a,b)
    #p0=function(x) {gamma(x[2]+x[4])*gamma(x[3]+x[4])/(gamma(x[2]+x[3]+x[4])*gamma(x[4]));}; #x[1]=n
    p02 = function(x){extraDistr::dbnbinom(0, size = x[2], alpha = x[3], beta = x[4]);}; #c(0, n, a,b)
    pthetad.r = function(x){extraDistr::dbnbinom(x[1], size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[1]) + digamma(x[2]+x[3]) - digamma(x[2]) - digamma(x[1]+x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbnbinom(x[1], size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) + digamma(x[3]+x[4]) - digamma(x[1]+x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbnbinom(x[1], size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[1]+x[3]+x[2]+x[4]) - digamma(x[4]));};
    p0d.r = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) + digamma(x[3]+x[4]) - digamma(x[2]+x[3]+x[4]) + digamma(x[3]));};
    p0d.b = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    if(is.null(b0)) {b0=rep(0.3,4*np);}
    neg.log.lik <- function(theta) {
      g  <- theta[1:np];              # gamma
      b1 <- theta[np+(1:np)];         # beta1
      b2 <- theta[2*np+(1:np)];       # beta2
      b3 <- theta[3*np+(1:np)];       # beta3
      xg = x%*%g;                     # X^T*gamma
      xb1 = x%*%b1;                   # X^T*beta1
      xb2 = x%*%b2;                   # X^T*beta2
      xb3 = x%*%b3;                   # X^T*beta3
      phii = ginv(xg);                # g-inverse
      theta1i = h1inv(xb1)            # theta_1 i
      theta2i = h2inv(xb2)            # theta_2 i
      theta3i = h3inv(xb3)            # theta_3 i
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #if((min(thetay[,2]-thetay[,1])>=0)&(min(thetai0[,2])>0)&(min(thetai0[,3])>0)&(min(thetai0[,4])>0)&(min(phii)>0)&(max(phii)<1)) {
      #if((min(thetai0[,2])>0)&(min(thetai0[,3])>0)&(min(thetai0[,4])>0)&(min(phii)>0)&(max(phii)<1)) {
      ans = -sum(log(phii)[y==0]) - sum(log(1-phii)[y>0]) -
        sum(log(apply(thetay1, 1, ptheta2))) + sum(log(1 - apply(thetai01, 1, p02)));
      #};
      return(ans)
    }
    gp <- function(theta) {
      g  <- theta[1:np];               # gamma
      b1 <- theta[np+(1:np)];          # beta1
      b2 <- theta[2*np+(1:np)];        # beta2
      b3 <- theta[3*np+(1:np)];        # beta3
      xg = x%*%g;                      # X^T*gamma
      xb1 = x%*%b1;                    # X^T*beta1
      xb2 = x%*%b2;                    # X^T*beta2
      xb3 = x%*%b3;                    # X^T*beta3
      phii = ginv(xg);                 # g-inverse
      theta1i =  ceiling(h1inv(xb1))   # h1-inverse
      theta2i = h2inv(xb2)             # h2-inverse
      theta3i = h3inv(xb3)             # h3-inverse
      phiid = ginvd(xg)                # g-inverse derivative
      theta1id = h1invd(xb1)           # h1-inverse derivative
      theta2id = h2invd(xb2)           # h2-inverse derivative
      theta3id = h3invd(xb3)           # h3-inverse derivative
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00=thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      thetai01 = thetai0[y>0,];
      #wrt phi
      dl = -apply((phiid/(phii))[y==0] * x[y==0,], 2, sum) + apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      # wrt r
      c1a = apply(thetay1, 1, pthetad.r) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta2)
      c1c = (c1a/c1b) * x[y>0,]
      temp3.r = apply(c1c, 2, sum)
      d1a = apply(thetai01, 1, p0d.r) * theta1id[y>0]
      d1b = 1 - (apply(thetai01, 1, p02))
      d1c = (d1a/d1b) * x[y>0,]
      temp4.r = apply(d1c, 2, sum)
      dr = -temp3.r - temp4.r
      #wrt a
      e1 = apply(thetay1, 1, pthetad.a) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta2)
      e3 = (e1/e2) * x[y>0,]
      temp3.a = apply(e3, 2, sum)
      f1 = apply(thetai01, 1, p0d.a) * theta2id[y>0]
      f2 = 1 - (apply(thetai01, 1, p02))
      f3 = (f1/f2) * x[y>0,]
      temp4.a = apply(f3, 2, sum)
      da = -temp3.a - temp4.a
      #wrt b
      z1 = apply(thetay1, 1, pthetad.b) * (theta3id[y>0])
      z2 = apply(thetay1, 1, ptheta2)
      z3 = (z1/z2) * x[y>0,]
      temp3.b = apply(z3, 2, sum)
      k1 = apply(thetai01, 1, p0d.b) * theta3id[y>0]
      k2 = 1 - (apply(thetai01, 1, p02))
      k3 = (k1/k2) * x[y>0,]
      temp4.b = apply(k3, 2, sum)
      db = -temp3.b - temp4.b
      return(c(dl, dr, da, db))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method =  "BFGS")
    estlike = estimate$value
    Akaike.information.criterion = -2 * -(estlike) + 2 * 4 * np
    Bayesian.information.criterion = -2 * -(estlike) + (log (N)) * 4 * np
    mle = matrix(c(estimate$par), nrow = 1)
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta1 <- mle[vec2]
    vec3 <- 2*np+(1:np)
    beta2 <- mle[vec3]
    vec4 <- 3*np+(1:np)
    beta3 <- mle[vec4]
    Estimated = matrix (c(gamma, beta1, beta2, beta3), ncol = np, nrow = 4, byrow = TRUE)
    rownames(Estimated) = c("Estimated Gamma", "Estimated Beta1", "Estimated Beta2", "Estimated Beta3")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike , Estimated.Parameters = Estimated ))
  }
  if (dist == "ZIBB")    {
    h1inv = function(theta) {exp(theta);};    # h1-inverse
    h2inv = function(theta) {exp(theta);};    # h2-inverse
    h3inv = function(theta) {exp(theta);};    # h3-inverse
    h1invd = function(theta) {exp(theta);};   # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};   # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};   # h3-inverse derivative
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};   # c(y, n, a, b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    neg.log.lik <- function(theta) {
      g  <- theta[1:np];              # gamma
      b1 <- theta[np+(1:np)];         # beta1
      b2 <- theta[2*np+(1:np)];       # beta2
      b3 <- theta[3*np+(1:np)];       # beta3
      xg = x%*%g;                     # X^T*gamma
      xb1 = x%*%b1;                   # X^T*beta1
      xb2 = x%*%b2;                   # X^T*beta2
      xb3 = x%*%b3;                   # X^T*beta3
      phii = ginv(xg);                # g-inverse
      theta1i = ceiling(h1inv(xb1))   # theta_1 i
      theta2i = h2inv(xb2)            # theta_2 i
      theta3i = h3inv(xb3)            # theta_3 i
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      ans = -sum(log(phii[y == 0] + (1 - phii[y == 0]) * apply(thetai00, 1, p03))) -
             sum(log(1 - phii)[y > 0]) - sum(log(apply(thetay1, 1, ptheta3)));
      return(ans)
    }
    gp <- function(theta) {
      g  <- theta[1:np];               # gamma
      b1 <- theta[np+(1:np)];          # beta1
      b2 <- theta[2*np+(1:np)];        # beta2
      b3 <- theta[3*np+(1:np)];        # beta3
      xg = x%*%g;                      # X^T*gamma
      xb1 = x%*%b1;                    # X^T*beta1
      xb2 = x%*%b2;                    # X^T*beta2
      xb3 = x%*%b3;                    # X^T*beta3
      phii = ginv(xg);                 # g-inverse
      theta1i =  ceiling(h1inv(xb1))   # h1-inverse
      theta2i = h2inv(xb2)             # h2-inverse
      theta3i = h3inv(xb3)             # h3-inverse
      phiid = ginvd(xg)                # g-inverse derivative
      theta1id = h1invd(xb1)           # h1-inverse derivative
      theta2id = h2invd(xb2)           # h2-inverse derivative
      theta3id = h3invd(xb3)           # h3-inverse derivative
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #wrt phi
      #temp1 = apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      temp1 <- colSums((phiid / (1 - phii))[y > 0] * x[y > 0, ])
      temp2a = (1 - apply(thetai00, 1, p03)) * phiid[y==0]
      temp2b = phii[y==0] + ((1 - phii)[y==0]) * apply(thetai00, 1, p03)
      temp2c <- (temp2a/temp2b) * x[y==0,]
      #temp2 = apply(temp2c, 2, sum)
      temp2 <- colSums(temp2c)
      dl = temp1 - temp2
      #wrt n
      c1a = apply(thetay1, 1, pthetad.n) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta3)
      c1c = (c1a/c1b) * x[y>0,]
      temp3.n = apply(c1c, 2, sum)
      d1a = (1 - phii[y==0]) * apply(thetai00, 1, p0d.n) * theta1id[y==0]
      d1b = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      d1c = (d1a/d1b) * x[y==0,]
      temp4.n = apply(d1c, 2, sum)
      dn = -temp3.n - temp4.n
      #wrt a
      e1 = apply(thetay1, 1, pthetad.a) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta3)
      e3 = (e1/e2) * x[y>0,]
      temp3.a = apply(e3, 2, sum)
      f1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.a) * theta2id[y==0]
      f2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      f3 = (f1/f2) * x[y==0,]
      temp4.a = apply(f3, 2, sum)
      da = -temp3.a - temp4.a
      #wrt b
      z1 = apply(thetay1, 1, pthetad.b) * (theta3id[y>0])
      z2 = apply(thetay1, 1, ptheta3)
      z3 = (z1/z2) * x[y>0,]
      temp3.b = apply(z3, 2, sum)
      k1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.b) * theta3id[y==0]
      k2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      k3 = (k1/k2) * x[y==0,]
      temp4.b = apply(k3, 2, sum)
      db = -temp3.b - temp4.b
      #################################
      return(c(dl, dn, da, db))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    #estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "L-BFGS-B",lower = c(-Inf, rep(0, 2*np), rep(0, 2*np)), upper = rep(Inf, 4*np))
    estlike = estimate$value
    mle = matrix(c(estimate$par), nrow = 1)
    Akaike.information.criterion = -2 * -(estlike) + 2 * 4 * np
    Bayesian.information.criterion = -2 * -(estlike) + (log (N)) * 4 * np
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta1 <- mle[vec2]
    vec3 <- 2*np+(1:np)
    beta2 <- mle[vec3]
    vec4 <- 3*np+(1:np)
    beta3 <- mle[vec4]
    Estimated = matrix (c(gamma, beta1, beta2, beta3), ncol = np, nrow = 4, byrow = TRUE)
    rownames(Estimated) = c("Estimated Gamma", "Estimated Beta1", "Estimated Beta2", "Estimated Beta3")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike , Estimated.Parameters = Estimated ))
  }
  if (dist == "ZIBB-ab") {
    h1inv = function(theta) {exp(theta);};    # h1-inverse
    h2inv = function(theta) {exp(theta);};    # h2-inverse
    h3inv = function(theta) {exp(theta);};    # h3-inverse
    h1invd = function(theta) {exp(theta);};   # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};   # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};   # h3-inverse derivative
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};   # c(y, n, a, b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N))
    neg.log.lik <- function(theta) {
      g <- theta[1:np];        # gamma
      b1 <- theta[np+(1:np)];  # beta1
      b2 <- theta[2*np+1];     # beta2
      b3 <- theta[2*np+2];     # beta3
      xg = x%*%g;              # X^T*gamma
      xb1 = x%*%b1;            # X^T*beta1
      xb2 = x1%*%b2;           # X^T*beta2
      xb3 = x1%*%b3;           # X^T*beta3
      phii = ginv(xg);         # phi_i
      theta1i = ceiling(h1inv(xb1))
      theta2i = h2inv(xb2)
      theta3i = h3inv(xb3)
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
        ans = -sum(log(1 - phii)[y>0]) - sum(log(apply(thetay1, 1, ptheta3))) -
               sum(log(phii[y==0] + (1 - phii[y==0]) * apply(thetai00, 1, p03)));
    }
    gp <- function(theta) {
      g <- theta[1:np];        # gamma
      b1 <- theta[np+(1:np)];  # beta1
      b2 <- theta[2*np+1];     # beta2
      b3 <- theta[2*np+2];     # beta3
      xg = x%*%g;              # X^T*gamma
      xb1 = x%*%b1;            # X^T*beta1
      xb2 = x1%*%b2;           # X^T*beta2
      xb3 = x1%*%b3;           # X^T*beta3
      phii = ginv(xg);         # phi_i
      theta1i = ceiling(h1inv(xb1))
      theta2i = h2inv(xb2)
      theta3i = h3inv(xb3)
      phiid = ginvd(xg)
      theta1id = h1invd(xb1)
      theta2id = h2invd(xb2)
      theta3id = h3invd(xb3)
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00=thetai0[y==0,];
      # wrt phi
      temp1 = apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      temp2a = (1 - apply(thetai00, 1, p03)) * phiid[y==0]
      temp2b = phii[y==0] + ((1 - phii)[y==0]) * apply(thetai00, 1, p03)
      temp2c = (temp2a/temp2b) * x[y==0,]
      temp2 = apply(temp2c, 2, sum)
      dl = temp1 - temp2
      # wrt n
      c1a = apply(thetay1, 1, pthetad.n) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta3)
      c1c = (c1a/c1b) * x[y>0,]
      temp3.n = apply(c1c, 2, sum)
      d1a = (1 - phii[y==0]) * apply(thetai00, 1, p0d.n) * theta1id[y==0]
      d1b = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      d1c = (d1a/d1b) * x[y==0,]
      temp4.n = apply(d1c, 2, sum)
      dn = -temp3.n - temp4.n
      # wrt a
      e1 = apply(thetay1, 1, pthetad.a) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta3)
      e3 = as.matrix(((e1/e2)*x1[y>0,]), ncol=4)
      temp3.a=sum(e3)
      f1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.a) * theta2id[y==0]
      f2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      f3 = (f1/f2) * x1[y==0,]
      temp4.a = sum(f3)
      da = -temp3.a - temp4.a
      # wrt b
      z1 = apply(thetay1, 1, pthetad.b) * (theta3id[y>0])
      z2  =apply(thetay1, 1, ptheta3)
      z3 = (z1/z2) * x1[y>0,]
      temp3.b = sum(z3)
      k1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.b) * theta3id[y==0]
      k2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      k3 = (k1/k2) * x1[y==0,]
      temp4.b = sum(k3)
      db = -temp3.b - temp4.b
      #################################
      return(c(dl, dn, da, db))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    Akaike.information.criterion = -2 * -(estlike) + 2 * 4 * np
    Bayesian.information.criterion = -2 * -(estlike) + (log (N)) * 4 * np
    mle = matrix(c(estimate$par), nrow = 1)
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta1 <- mle[vec2]
    vec3 <- 2*np+1
    beta2 <- mle[vec3]
    vec4 <- 2*np+2
    beta3 <- mle[vec4]
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike ,
                Estimated.Gamma = gamma, Estimated.Beta1 = beta1, Estimated.Beta2 = beta2, Estimated.Beta3 = beta3 ))
  }
  if (dist == "ZIBB-n")  {
    h1inv = function(theta) {exp(theta);};    # h1-inverse
    h2inv = function(theta) {exp(theta);};    # h2-inverse
    h3inv = function(theta) {exp(theta);};    # h3-inverse
    h1invd = function(theta) {exp(theta);};   # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};   # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};   # h3-inverse derivative
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};   # c(y, n, a, b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N))
    neg.log.lik <- function(theta) {
      g <- theta[1:np];           # gamma
      b1 <- theta[np+1];          # beta1
      b2 <- theta[np+1+(1:np)];   # beta2
      b3 <- theta[2*np+1+(1:np)]; # beta3
      xg = x%*%g;                 # X^T*gamma
      xb1 = x1%*%b1;              # X^T*beta1
      xb2 = x%*%b2;               # X^T*beta2
      xb3 = x%*%b3;               # X^T*beta3
      phii = ginv(xg);            #phi_i
      theta1i = ceiling(h1inv(xb1))
      theta2i = h2inv(xb2)
      theta3i = h3inv(xb3)
      thetai0=cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00=thetai0[y==0,];
      ans = -sum(log(1 - phii)[y>0]) - sum(log(apply(thetay1, 1, ptheta3))) -
             sum(log(phii[y==0] + (1 - phii[y==0]) * apply(thetai00, 1, p03)));
      return(ans)
    }
    gp <- function(theta) {
      g <- theta[1:np];                 # gamma
      b1 <- theta[np+1];                # beta1
      b2 <- theta[np+1+(1:np)];         # beta2
      b3 <- theta[2*np+1+(1:np)];       # beta3
      xg = x%*%g;                       # X^T*gamma
      xb1 = x1%*%b1;                    # X^T*beta1
      xb2 = x%*%b2;                     # X^T*beta2
      xb3 = x%*%b3;                     # X^T*beta3
      phii = ginv(xg);                  # g-inverse
      theta1i =  ceiling(h1inv(xb1))    # h1-inverse
      theta2i = h2inv(xb2)              # h2-inverse
      theta3i = h3inv(xb3)              # h3-inverse
      phiid = ginvd(xg)                 # g-inverse derivative
      theta1id = h1invd(xb1)            # h1-inverse derivative
      theta2id = h2invd(xb2)            # h2-inverse derivative
      theta3id = h3invd(xb3)            # h3-inverse derivative
      thetai0=cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00=thetai0[y==0,];
      # wrt phi
      temp1 = apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      temp2a = (1 - apply(thetai00, 1, p03)) * phiid[y==0]
      temp2b = phii[y==0] + ((1 - phii)[y==0]) * apply(thetai00, 1, p03)
      temp2c = (temp2a/temp2b) * x[y==0,]
      temp2 = apply(temp2c, 2, sum)
      dl = temp1 - temp2
      # wrt n
      c1a = apply(thetay1, 1, pthetad.n) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta3)
      c1c = (c1a/c1b )* x1[y>0,]
      temp3.n = sum(c1c)
      d1a = (1 - phii[y==0]) * apply(thetai00, 1, p0d.n) * theta1id[y==0]
      d1b = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      d1c = (d1a/d1b) * x1[y==0,]
      temp4.n = sum(d1c)
      dn = -temp3.n - temp4.n
      # wrt a
      e1 = apply(thetay1, 1, pthetad.a) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta3)
      e3 = as.matrix(((e1/e2)*x[y>0,]), ncol=4)
      temp3.a = apply(e3, 2, sum)
      f1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.a) * theta2id[y==0]
      f2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      f3 = (f1/f2) * x[y==0,]
      temp4.a = apply(f3, 2, sum)
      da = -temp3.a - temp4.a
      # wrt b
      z1 = apply(thetay1, 1, pthetad.b) * (theta3id[y>0])
      z2 = apply(thetay1, 1, ptheta3)
      z3 = (z1/z2) * x[y>0,]
      temp3.b = apply(z3, 2, sum)
      k1 = (1 - phii[y==0]) * apply(thetai00, 1, p0d.b) * theta3id[y==0]
      k2 = phii[y==0] + (1 - phii[y==0]) * (apply(thetai00, 1, p03))
      k3 = (k1/k2) * x[y==0,]
      temp4.b = apply(k3, 2, sum)
      db = -temp3.b - temp4.b
      return(c(dl, dn, da, db))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    Akaike.information.criterion = -2 * -(estlike) + 2 * 4 * np
    Bayesian.information.criterion = -2 * -(estlike) + (log (N)) * 4 * np
    mle = matrix(c(estimate$par), nrow = 1)
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+1
    beta1 <- mle[vec2]
    vec3 <- np+1+(1:np)
    beta2 <- mle[vec3]
    vec4 <- 2*np+1+(1:np)
    beta3 <- mle[vec4]
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike ,
                Estimated.Gamma = gamma, Estimated.Beta1 = beta1, Estimated.Beta2 = beta2, Estimated.Beta3 = beta3 ))
  }
  if (dist == "BBH")     {
    h1inv = function(theta) {exp(theta);};    # h1-inverse
    h2inv = function(theta) {exp(theta);};    # h2-inverse
    h3inv = function(theta) {exp(theta);};    # h3-inverse
    h1invd = function(theta) {exp(theta);};   # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};   # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};   # h3-inverse derivative
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};   # c(y, n, a, b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    if(is.null(b0))  {b0 = rep(0.2, 4*np);}
    neg.log.lik <- function(theta) {
      g  <- theta[1:np];              # gamma
      b1 <- theta[np+(1:np)];         # beta1
      b2 <- theta[2*np+(1:np)];       # beta2
      b3 <- theta[3*np+(1:np)];       # beta3
      xg = x %*% g;                   # X^T*gamma
      xb1 = x %*% b1;                 # X^T*beta1
      xb2 = x %*% b2;                 # X^T*beta2
      xb3 = x %*% b3;                 # X^T*beta3
      phii = ginv(xg);                # g-inverse
      theta1i = ceiling(h1inv(xb1))   # theta_1 i
      #theta1i = h1inv(xb1)   # theta_1 i
      theta2i = h2inv(xb2)            # theta_2 i
      theta3i = h3inv(xb3)            # theta_3 i
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      #thetai0 = cbind(theta1i, theta2i, theta3i);
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #if((min(thetay[,2]-thetay[,1])>=0)&(min(thetai0[,2])>0)&(min(thetai0[,3])>0)&(min(thetai0[,4])>0)&(min(phii)>0)&(max(phii)<1)) {
      # if((min(thetai0[,2])>0)&(min(thetai0[,3])>0)&(min(thetai0[,4])>0)&(min(phii)>0)&(max(phii)<1)) {
      ans = -sum(log(phii)[y==0]) - sum(log(1 - phii)[y>0]) -
        sum(log(apply(thetay1, 1, ptheta3))) + sum(log(1 - apply(thetai01, 1, p03)));
      #};
      return(ans)
    }
    gp <- function(theta) {
      g  <- theta[1:np];               # gamma
      b1 <- theta[np+(1:np)];          # beta1
      b2 <- theta[2*np+(1:np)];        # beta2
      b3 <- theta[3*np+(1:np)];        # beta3
      xg = x%*%g;                      # X^T*gamma
      xb1 = x%*%b1;                    # X^T*beta1
      xb2 = x%*%b2;                    # X^T*beta2
      xb3 = x%*%b3;                    # X^T*beta3
      phii = ginv(xg);                 # g-inverse
      theta1i =  ceiling(h1inv(xb1))   # h1-inverse
      #theta1i =  h1inv(xb1)   # h1-inverse
      theta2i = h2inv(xb2)             # h2-inverse
      theta3i = h3inv(xb3)             # h3-inverse
      phiid = ginvd(xg)                # g-inverse derivative
      theta1id = h1invd(xb1)           # h1-inverse derivative
      theta2id = h2invd(xb2)           # h2-inverse derivative
      theta3id = h3invd(xb3)           # h3-inverse derivative
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      #thetai0 = cbind(theta1i, theta2i, theta3i);
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #wrt phi
      dl = -apply((phiid/(phii))[y==0] * x[y==0,], 2, sum) + apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      #wrt n
      c1a = apply(thetay1, 1, pthetad.n) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta3)
      c1c = (c1a/c1b) * x[y>0,]
      temp3.n = apply(c1c, 2, sum)
      d1a = apply(thetai01, 1, p0d.n) * theta1id[y>0]
      d1b = 1 - (apply(thetai01, 1, p03))
      d1c = (d1a/d1b) * x[y>0,]
      temp4.n = apply(d1c, 2, sum)
      dn = -temp3.n - temp4.n
      #wrt a
      e1 = apply(thetay1, 1, pthetad.a) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta3)
      e3 = (e1/e2) * x[y>0,]
      temp3.a = apply(e3, 2, sum)
      f1 = apply(thetai01, 1, p0d.a) * theta2id[y>0]
      f2 = 1 - (apply(thetai01, 1, p03))
      f3 = (f1/f2) * x[y>0,]
      temp4.a = apply(f3, 2, sum)
      da = -temp3.a - temp4.a
      #wrt b
      z1 = apply(thetay1, 1, pthetad.b) * (theta3id[y>0])
      z2 = apply(thetay1, 1, ptheta3)
      z3 = (z1/z2) * x[y>0,]
      temp3.b = apply(z3, 2, sum)
      k1 = apply(thetai01, 1, p0d.b) * theta3id[y>0]
      k2 = 1 - (apply(thetai01, 1, p03))
      k3 = (k1/k2) * x[y>0,]
      temp4.b = apply(k3, 2, sum)
      db = -temp3.b - temp4.b
      #################################
      return(c(dl, dn, da, db))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    Akaike.information.criterion = -2 * -(estlike) + 2 * 4 * np
    Bayesian.information.criterion = -2 * -(estlike) + log (N) * 4 * np
    mle = matrix(c(estimate$par), nrow = 1)
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta1 <- mle[vec2]
    vec3 <- 2*np+(1:np)
    beta2 <- mle[vec3]
    vec4 <- 3*np+(1:np)
    beta3 <- mle[vec4]
    Estimated = matrix (c(gamma, beta1, beta2, beta3), ncol = np, nrow = 4, byrow = TRUE)
    rownames(Estimated) = c("Estimated Gamma", "Estimated Beta1", "Estimated Beta2", "Estimated Beta3")
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike , Estimated.Parameters = Estimated ))
  }
  if (dist == "BBH-ab")  {
    h1inv = function(theta) {exp(theta);};    # h1-inverse
    h2inv = function(theta) {exp(theta);};    # h2-inverse
    h3inv = function(theta) {exp(theta);};    # h3-inverse
    h1invd = function(theta) {exp(theta);};   # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};   # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};   # h3-inverse derivative
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};   # c(y, n, a, b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N))
    #if(is.null(b0))  {b0=rep(0.2,4*np);}
    neg.log.lik <- function(theta) {
      g  <- theta[1:np];              # gamma
      b1 <- theta[np+(1:np)];         # beta1
      b2 <- theta[2*np+1];            # beta2
      b3 <- theta[2*np+2];            # beta3
      xg = x%*%g;                     # X^T*gamma
      xb1 = x%*%b1;                   # X^T*beta1
      xb2 = x1%*%b2;                  # X^T*beta2
      xb3 = x1%*%b3;                  # X^T*beta3
      phii = ginv(xg);                # g-inverse
      theta1i = ceiling(h1inv(xb1))   # theta_1 i
      theta2i = h2inv(xb2)            # theta_2 i
      theta3i = h3inv(xb3)            # theta_3 i
      thetai0=cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #loglik=-Inf;
      #if((min(thetay[,2]-thetay[,1])>=0)&(min(thetai0[,2])>0)&(min(thetai0[,3])>0)&(min(thetai0[,4])>0)&(min(phii)>0)&(max(phii)<1)) {
      # if((min(thetai0[,2])>0)&(min(thetai0[,3])>0)&(min(thetai0[,4])>0)&(min(phii)>0)&(max(phii)<1)) {
      ans = -sum(log(phii)[y==0]) - sum(log(1-phii)[y>0]) -
        sum(log(apply(thetay1, 1, ptheta3))) + sum(log(1 - apply(thetai01, 1, p03)));
      #};
      return(ans)
    }
    gp <- function(theta) {
      g  <- theta[1:np];               # gamma
      b1 <- theta[np+(1:np)];          # beta1
      b2 <- theta[2*np+1];             # beta2
      b3 <- theta[2*np+2];             # beta3
      xg = x%*%g;                      # X^T*gamma
      xb1 = x%*%b1;                    # X^T*beta1
      xb2 = x1%*%b2;                   # X^T*beta2
      xb3 = x1%*%b3;                   # X^T*beta3
      phii = ginv(xg);                 # g-inverse
      theta1i =  ceiling(h1inv(xb1))   # h1-inverse
      theta2i = h2inv(xb2)             # h2-inverse
      theta3i = h3inv(xb3)             # h3-inverse
      phiid = ginvd(xg)                # g-inverse derivative
      theta1id = h1invd(xb1)           # h1-inverse derivative
      theta2id = h2invd(xb2)           # h2-inverse derivative
      theta3id = h3invd(xb3)           # h3-inverse derivative
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #wrt phi
      dl = -apply((phiid/(phii))[y==0] * x[y==0,], 2, sum) + apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      #wrt n
      c1a = apply(thetay1, 1, pthetad.n) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta3)
      c1c = (c1a/c1b) * x[y>0,]
      temp3.n = apply(c1c, 2, sum)
      d1a = apply(thetai01, 1, p0d.n) * theta1id[y>0]
      d1b = 1 - (apply(thetai01, 1, p03))
      d1c = (d1a/d1b) * x[y>0,]
      temp4.n = apply(d1c, 2, sum)
      dn = -temp3.n - temp4.n
      #wrt a
      e1 = apply(thetay1, 1, pthetad.a) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta3)
      e3 = as.matrix(((e1/e2) * x1[y>0,]), ncol=4)
      temp3.a = sum(e3)
      f1 = apply(thetai01, 1, p0d.a) * theta2id[y>0]
      f2 = 1 - (apply(thetai01, 1, p03))
      f3 = (f1/f2) * x1[y>0,]
      temp4.a = sum(f3)
      da = -temp3.a - temp4.a
      #wrt b
      z1 = apply(thetay1, 1, pthetad.b) * (theta3id[y>0])
      z2 = apply(thetay1, 1, ptheta3)
      z3 = (z1/z2) * x1[y>0,]
      temp3.b = sum(z3)
      k1 = apply(thetai01, 1, p0d.b) * theta3id[y>0]
      k2 = 1 - (apply(thetai01, 1, p03))
      k3 = (k1/k2) * x1[y>0,]
      temp4.b = sum(k3)
      db = -temp3.b - temp4.b
      #################################
      return(c(dl, dn, da, db))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    Akaike.information.criterion = -2 * -(estlike) + 2 * 4 * np
    Bayesian.information.criterion = -2 * -(estlike) + log (N) * 4 * np
    mle = matrix(c(estimate$par), nrow = 1)
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+(1:np)
    beta1 <- mle[vec2]
    vec3 <- 2*np+1
    beta2 <- mle[vec3]
    vec4 <- 2*np+2
    beta3 <- mle[vec4]
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike ,
                Estimated.Gamma = gamma, Estimated.Beta1 = beta1, Estimated.Beta2 = beta2, Estimated.Beta3 = beta3 ))
  }
  if (dist == "BBH-n")   {
    h1inv = function(theta) {exp(theta);};    # h1-inverse
    h2inv = function(theta) {exp(theta);};    # h2-inverse
    h3inv = function(theta) {exp(theta);};    # h3-inverse
    h1invd = function(theta) {exp(theta);};   # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};   # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};   # h3-inverse derivative
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};   # c(y, n, a, b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    x1 = as.matrix(rep(1, N))
    #if(is.null(b0))  {b0=rep(0.2,4*np);}
    neg.log.lik <- function(theta) {
      g  <- theta[1:np];              # gamma
      b1 <- theta[np+1];              # beta1
      b2 <- theta[np+1+(1:np)];       # beta2
      b3 <- theta[2*np+1+(1:np)];     # beta3
      xg = x%*%g;                     # X^T*gamma
      xb1 = x1%*%b1;                  # X^T*beta1
      xb2 = x%*%b2;                   # X^T*beta2
      xb3 = x%*%b3;                   # X^T*beta3
      phii = ginv(xg);                # g-inverse
      theta1i = ceiling(h1inv(xb1))   # theta_1 i
      theta2i = h2inv(xb2)            # theta_2 i
      theta3i = h3inv(xb3)            # theta_3 i
      thetai0=cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #loglik=-Inf;
      #if((min(thetay[,2]-thetay[,1])>=0)&(min(thetai0[,2])>0)&(min(thetai0[,3])>0)&(min(thetai0[,4])>0)&(min(phii)>0)&(max(phii)<1)) {
      # if((min(thetai0[,2])>0)&(min(thetai0[,3])>0)&(min(thetai0[,4])>0)&(min(phii)>0)&(max(phii)<1)) {
      ans = -sum(log(phii)[y==0]) - sum(log(1-phii)[y>0]) -
        sum(log(apply(thetay1, 1, ptheta3))) + sum(log(1 - apply(thetai01, 1, p03)));
      #};
      return(ans)
    }
    gp <- function(theta) {
      g  <- theta[1:np];               # gamma
      b1 <- theta[np+1];               # beta1
      b2 <- theta[np+1+(1:np)];        # beta2
      b3 <- theta[2*np+1+(1:np)];      # beta3
      xg = x%*%g;                      # X^T*gamma
      xb1 = x1%*%b1;                   # X^T*beta1
      xb2 = x%*%b2;                    # X^T*beta2
      xb3 = x%*%b3;                    # X^T*beta3
      phii = ginv(xg);                 # g-inverse
      theta1i =  ceiling(h1inv(xb1))   # h1-inverse
      theta2i = h2inv(xb2)             # h2-inverse
      theta3i = h3inv(xb3)             # h3-inverse
      phiid = ginvd(xg)                # g-inverse derivative
      theta1id = h1invd(xb1)           # h1-inverse derivative
      theta2id = h2invd(xb2)           # h2-inverse derivative
      theta3id = h3invd(xb3)           # h3-inverse derivative
      thetai0 = cbind(0, theta1i, theta2i, theta3i)
      thetay = cbind(y, theta1i, theta2i, theta3i);
      thetay1 = thetay[y>0,];
      thetay0 = thetay[y==0,];
      thetai00 = thetai0[y==0,];
      thetai01 = thetai0[y>0,];
      #wrt phi
      dl = -apply((phiid/(phii))[y==0] * x[y==0,], 2, sum) + apply((phiid/(1 - phii))[y>0] * x[y>0,], 2, sum)
      #wrt n
      c1a = apply(thetay1, 1, pthetad.n) * (theta1id[y>0])
      c1b = apply(thetay1, 1, ptheta3)
      c1c = (c1a/c1b) * x1[y>0,]
      temp3.n = sum(c1c)
      d1a = apply(thetai01, 1, p0d.n) * theta1id[y>0]
      d1b = 1 - (apply(thetai01, 1, p03))
      d1c = (d1a/d1b) * x1[y>0,]
      temp4.n = sum(d1c)
      dn = -temp3.n - temp4.n
      #wrt a
      e1 = apply(thetay1, 1, pthetad.a) * (theta2id[y>0])
      e2 = apply(thetay1, 1, ptheta3)
      e3 = (e1/e2) * x[y>0,]
      temp3.a = apply(e3, 2, sum)
      f1 = apply(thetai01, 1, p0d.a) * theta2id[y>0]
      f2 = 1 - (apply(thetai01, 1, p03))
      f3 = (f1/f2) * x[y>0,]
      temp4.a = apply(f3, 2, sum)
      da = -temp3.a - temp4.a
      #wrt b
      z1 = apply(thetay1, 1, pthetad.b) * (theta3id[y>0])
      z2 = apply(thetay1, 1, ptheta3)
      z3 = (z1/z2) * x[y>0,]
      temp3.b = apply(z3, 2, sum)
      k1 = apply(thetai01, 1, p0d.b) * theta3id[y>0]
      k2 = 1 - (apply(thetai01, 1, p03))
      k3 = (k1/k2) * x[y>0,]
      temp4.b = apply(k3, 2, sum)
      db = -temp3.b - temp4.b
      #################################
      return(c(dl, dn, da, db))
    }
    estimate = stats::optim(par = b0, fn = neg.log.lik, gr = gp, method = "BFGS")
    estlike = estimate$value
    Akaike.information.criterion = -2 * -(estlike) + 2 * 4 * np
    Bayesian.information.criterion = -2 * -(estlike) + log (N) * 4 * np
    mle = matrix(c(estimate$par), nrow = 1)
    vec1 <- 1:np
    gamma <- mle[vec1]
    vec2 <- np+1
    beta1 <- mle[vec2]
    vec3 <- np+1+(1:np)
    beta2 <- mle[vec3]
    vec4 <- 2*np+1+(1:np)
    beta3 <- mle[vec4]
    return(list(AIC = Akaike.information.criterion, BIC = Bayesian.information.criterion, log.likelihhod = estlike ,
                Estimated.Gamma = gamma, Estimated.Beta1 = beta1, Estimated.Beta2 = beta2, Estimated.Beta3 = beta3 ))
  }
}
