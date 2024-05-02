#' Fisher Information computation for calculating the confidence intervals in Zero-inflated and Zero-altered regression models
#'
#' @param x a design matrix containing an intercept column (all ones) along with other available covariates for the response variable
#' @param y a zero-inflated or zero-altered(Hurdle) count response variable, represented as an integer vector
#' @param b0 the initial parameters for the model, calculated as the product of the number of parameters in the specified models and the number of covariates
#' @param m M set in trigamma free approach only needed for ZIBNB, BNBH, ZINB, and NBH
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
#' @return
#' If dist = ZIP, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = ZINB, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = ZINB-r, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = ZIBNB, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = ZIBB, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = ZIBB-n, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = ZIBB-ab, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = PH, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = NBH, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = NBH-r, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = BNBH, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = BBH, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = BBH-n, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' If dist = BBH-ab, the following values are returned:
#' \itemize{
#' \item FisherInformation: Fisher Information matrix for all the parameters and covariates in the model.
#' \item ConfidenceIntervals: Contains the following information:
#'     \itemize{
#'       \item Lower and upper bounds of the confidence interval.
#'       \item Estimated parameters.
#'       \item Estimation/length ratio.
#'       \item Standard error.
#'       \item Z-score.
#'     }
#' }
#' @export
#' @examples
#'intercept<- rep(1,4406)
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
#' b0 = c (0.1,0.10,0.1,0.1,0.1,0.1,0.1,0.1,
#'0.2,0.20,0.20,0.20,0.20,0.20,0.20,0.20)
#'fisher.reg(x, y, b0=b0, m=1e3, dist="ZIP", link="probit")
fisher.reg <- function (x, y, b0 = NULL, m, dist = "ZIP", link = "logit"){
  if (link == "loglog")    {
    ginv = function(phi) {exp(-exp(phi));};                       # g-inverse with "log-log" link
    ginvd = function(phi){-exp(phi-exp(phi));};                   # g-inverse derivative
  }
  if (link == "cloglog")   {
    ginv = function(phi) {1-exp(-exp(phi));};                     # g-inverse with "c-log-log"
    ginvd = function(phi){exp(phi-exp(phi));};                    # g-inverse derivative
  }
  if (link == "logit")     {
    ginv = function(phi){exp(phi)/(1 + exp(phi));};               # g-inverse with "logit" link
    ginvd = function(phi){exp(phi)/(1+exp(phi))^2;};              # g-inverse derivative
  }
  if (link == "probit")    {
    ginv = function(phi){stats::pnorm(phi, mean = 0, sd = 1);};   # g-inverse with "probit" link
    ginvd = function(phi){stats::dnorm(phi, mean = 0, sd = 1);};  # g-inverse derivative
  }
  if (dist == "ZIP")       {
    ZIP <- RAZIAD::reg.model(x, y, b0=b0, dist="ZIP", link)
    g <- ZIP$Estimated.Parameters[1,]                            # Estimated gamma
    b <- ZIP$Estimated.Parameters[2,]                            # Estimated beta
    hinv = function(theta) {exp(theta);};                        # h-inverse
    hinvd = function(theta) {exp(theta);};                       # h-inverse derivative
    ptheta = function(x){stats::dpois(x[2], x[1]);};                    # c(theta, y)
    pthetad = function(x){stats::dpois(x[2], x[1])*(x[2]/x[1]-1);};     # c(theta,y)
    p0 = function(theta){exp(-theta);};
    p0d = function(theta){-exp(-theta);};
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    xg = x%*%g;                     # X^T*gamma
    xb = x%*%b;                     # X^T*beta1
    phii = ginv(xg);                # phi_i
    thetai = hinv(xb);              # theta_i
    p00 = exp(-hinv(xb))

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    C = as.vector(hinvd(xb)) * x
    #A11
    gg = t(A) %*% (A * as.vector((1/P) - 1))
    #A12
    gb0 = (1-phii) * (p00/P) * (-1)
    gb = t(A) %*% (C * as.vector(gb0))
    #A22
    eist = 1/(hinv(xb)) #1/lambda
    bb1 = t(C) %*% (C * as.vector((1-phii) * eist))
    H = (1/P) * phii * (1-phii) * p00
    bb2 = t(C) %*% (C * as.vector(H))
    bb = bb1 - bb2

    F1 = cbind(gg, gb)
    F2 = cbind(t(gb), bb )
    Fisher = round(rbind(F1, F2),2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }
    EST <- c(g,b)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*2), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt <- round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))}
  if (dist == "PH")        {
    PH <- RAZIAD::reg.model(x, y, b0=b0, dist="PH", link)
    g <- PH$Estimated.Parameters[1,]                             # Estimated gamma
    b <- PH$Estimated.Parameters[2,]                             # Estimated beta
    hinv = function(theta) {exp(theta);};                        # h-inverse
    hinvd = function(theta) {exp(theta);};                       # h-inverse derivative
    ptheta=function(x){stats::dpois(x[2], x[1]);};                      # c(theta, y)
    pthetad=function(x){stats::dpois(x[2], x[1])*(x[2]/x[1]-1);};
    p0=function(theta){exp(-theta);};
    p0d=function(theta){-exp(-theta);};
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    xg = x%*%g;                     # X^T*gamma
    xb = x%*%b;                     # X^T*beta1
    phii = ginv(xg);                # phi_i
    thetai = hinv(xb);              # theta_i

    P1 = 1/(phii * (1 - phii))
    P2 = (1 - phii)/(1 - p0(thetai))
    P3 = p0(thetai)/(1 - p0(thetai))
    A = as.vector(ginvd(xg)) * x
    C = as.vector(hinvd(xb)) * x

    #Fisher
    #A11
    gg = t(A)%*%(A * as.vector(P1))
    #A12
    gb = matrix(0, nrow = dim(x)[2], ncol = dim(x)[2])
    #A22
    eist = 1/(hinv(xb)) #1/lambda
    bb1 = t(C) %*% (C * as.vector((P2) * eist))
    bb2 = t(C) %*% (C * as.vector((1 - phii) * P3 * P2))
    bb = bb1 - bb2

    F1 = cbind(gg, gb)
    F2 = cbind(t(gb), bb)
    Fisher = round(rbind(F1, F2),2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }
    EST <- c(g,b)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*2), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt <- round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))}
  if (dist == "ZINB")      {
    ZINB <- RAZIAD::reg.model(x, y, b0=b0, dist="ZINB", link)
    g  <- ZINB$Estimated.Parameters[1,]                           # Estimated gamma
    b1 <- ZINB$Estimated.Parameters[2,]                           # Estimated beta1
    b2 <- ZINB$Estimated.Parameters[3,]                           # Estimated beta2
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size= y[2], prob=y[3]);};
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*(digamma(x[1]+x[2])+log(x[3])-digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*((x[2]/(x[3]))-(x[1]/(1-x[3])));};
    p01 = function(x) {x[2]^x[1];}; #x[2]=p
    p0d.r = function(x) {log(x[2])*((x[2])^x[1])}
    p0d.p = function(x) {x[1]*((x[2])^(x[1]-1))}
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    phii = ginv(xg);                # g-inverse
    theta1i = h1inv(xb1)            # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    thetai0 = cbind(0, theta1i, theta2i)
    p00 = apply(thetai0, 1, p01)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (log(theta2i))
    A12 = t(A) %*% (Ci1*as.vector(gb10))
    gb20 = (1-phii) * (p00/P) * (theta1i/theta2i)
    A13 = t(A) %*% (Ci2 * as.vector(gb20))
    #In this part we calculate B11
    ##Fisher Information Using the other definition and applying trigamma free Method
    #### self-defined function to calculate E Psi_1(nu + Y) given Y ~ NB(r, p)
    M = m
    fEPsi1BNB <- function(nu, r, p, M) {
      trigamma(nu) - sum(stats::pnbinom(0:M, size = r, prob = p, lower.tail = F)/(nu + (0:M))^2);
    }
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = E(y_i + r_i)
    for(i in 1:n){
      nu1 = theta1i[i]
      E11[i] = fEPsi1BNB(nu = nu1, r = theta1i[i], p = theta2i[i], M = M)
    }

    #B11: r*r
    b1b1.0 =  -(E11 + trigamma(theta1i)) * (1-phii)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (1/P) * phii * (1-phii) * p00 * (log(theta1i))^2
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #B12: r*p
    b1b2.0 = -(1/theta2i) * (1-phii)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M1 = (1/P) * phii * (1-phii) * p00 * ((log(theta1i)) * (theta1i/theta2i))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M1))
    B12 = b1b2.1 - M2
    #B22: p*p
    b1b4.0 = (theta1i/((theta2i)^2 * (1-theta2i))) * (1-phii)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (1/P) * phii * (1-phii) * p00 * ((theta1i/theta2i)^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2

    # Fisher Information Matrix
    F1 = cbind(  A11,    A12,    A13)
    F2 = cbind(t(A12),   B11,    B12)
    F3 = cbind(t(A13), t(B12),   B22)
    Fisher = round(rbind(F1, F2, F3), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }

    EST <- c(g,b1,b2)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "NBH")       {
    NBH <- RAZIAD::reg.model(x, y, b0=b0, dist="NBH", link)
    g  <- NBH$Estimated.Parameters[1,]                            # Estimated gamma
    b1 <- NBH$Estimated.Parameters[2,]                            # Estimated beta
    b2 <- NBH$Estimated.Parameters[3,]                            # Estimated beta
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size= y[2], prob=y[3]);};
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*(digamma(x[1]+x[2])+log(x[3])-digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*((x[2]/(x[3]))-(x[1]/(1-x[3])));};
    p01 = function(x) {x[2]^x[1];};
    p0d.r = function(x) {log(x[2])*((x[2])^x[1])}
    p0d.p = function(x) {x[1]*((x[2])^(x[1]-1))}
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    phii = ginv(xg);                # g-inverse
    theta1i = h1inv(xb1)            # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    thetai0 = cbind(0, theta1i, theta2i)
    p00 = apply(thetai0, 1, p01)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (log(theta2i))
    A12 = matrix(0, nrow = dim(x)[2], ncol = dim(x)[2])
    A13 = A12
    #In this part we calculate B11
    ##Fisher Information Usingthe other definition and applying trigamma free Method
    #### self-defined function to calculate E Psi_1(nu + Y) given Y ~ NB(r, p)
    M = m
    fEPsi2BNB <- function(nu, r, p, M) {
      trigamma(nu) - sum(stats::pnbinom(0:M, size = r, prob = p, lower.tail = F)/(nu + (0:M))^2);
    }
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = E(y_i + r_i)
    for(i in 1:n){
      nu1 = theta1i[i]
      E11[i] = fEPsi2BNB(nu = nu1, r = theta1i[i], p = theta2i[i], M = M)
    }
    #B11: r*r
    b1b1.0 =  -(E11 + trigamma(theta1i)) * ((1-phii)/(1-p00))
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (p00/(1-p00)) * (1-phii) * p00 * (log(theta1i))^2
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #B12: r*p
    b1b2.0 = -(1/theta2i) * ((1-phii)/(1-p00))
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M1 = (p00/(1-p00)) * (1-phii) * p00 * ((log(theta1i)) * (theta1i/theta2i))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M1))
    B12 = b1b2.1 - M2
    #B22: p*p
    b1b4.0 = (theta1i/((theta2i)^2 * (1-theta2i))) * ((1-phii)/(1-p00))
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (p00/(1-p00)) * (1-phii) * p00 * ((theta1i/theta2i)^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2

    # Fisher Information Matrix
    F1 = cbind(  A11,    A12,    A13)
    F2 = cbind(t(A12),   B11,    B12)
    F3 = cbind(t(A13), t(B12),   B22)
    Fisher = round(rbind(F1, F2, F3), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }

    EST <- c(g,b1,b2)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "ZINB-r")    {
    ZINB.r <- RAZIAD::reg.model(x, y, b0=b0, dist="ZINB-r", link)
    g  <- ZINB.r$Estimated.Parameters[1,]                           # Estimated gamma
    b1 <- ZINB.r$Estimated.Parameters[2,]                           # Estimated beta1
    b2 <- ZINB.r$Estimated.Parameters[3,]                           # Estimated beta2
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size= y[2], prob=y[3]);};
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*(digamma(x[1]+x[2])+log(x[3])-digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*((x[2]/(x[3]))-(x[1]/(1-x[3])));};
    p01 = function(x) {x[2]^x[1];}; #x[2]=p
    p0d.r = function(x) {log(x[2])*((x[2])^x[1])}
    p0d.p = function(x) {x[1]*((x[2])^(x[1]-1))}
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N));
    xg = x%*%g;                     # X^T*gamma
    xb1 = x1%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    phii = ginv(xg);                # g-inverse
    theta1i = h1inv(xb1)            # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    thetai0 = cbind(0, theta1i, theta2i)
    p00 = apply(thetai0, 1, p01)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x1
    Ci2 = as.vector(h2invd(xb2)) * x

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (log(theta2i))
    A12 = t(A) %*% (Ci1*as.vector(gb10))
    gb20 = (1-phii) * (p00/P) * (theta1i/theta2i)
    A13 = t(A) %*% (Ci2 * as.vector(gb20))
    #In this part we calculate B11
    ##Fisher Information Using the other definition and applying trigamma free Method
    #### self-defined function to calculate E Psi_1(nu + Y) given Y ~ NB(r, p)
    M = m
    fEPsi1BNB <- function(nu, r, p, M) {
      trigamma(nu) - sum(stats::pnbinom(0:M, size = r, prob = p, lower.tail = F)/(nu + (0:M))^2);
    }
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = E(y_i + r_i)
    for(i in 1:n){
      nu1 = theta1i[i]
      E11[i] = fEPsi1BNB(nu = nu1, r = theta1i[i], p = theta2i[i], M = M)
    }

    #B11: r*r
    b1b1.0 =  -(E11 + trigamma(theta1i)) * (1-phii)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (1/P) * phii * (1-phii) * p00 * (log(theta1i))^2
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #B12: r*p
    b1b2.0 = -(1/theta2i) * (1-phii)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M1 = (1/P) * phii * (1-phii) * p00 * ((log(theta1i)) * (theta1i/theta2i))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M1))
    B12 = b1b2.1 - M2
    #B22: p*p
    b1b4.0 = (theta1i/((theta2i)^2 * (1-theta2i))) * (1-phii)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (1/P) * phii * (1-phii) * p00 * ((theta1i/theta2i)^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2

    # Fisher Information Matrix
    F1 = cbind(  A11,    A12,    A13)
    F2 = cbind(t(A12),   B11,    B12)
    F3 = cbind(t(A13), t(B12),   B22)
    Fisher = round(rbind(F1, F2, F3), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }

    EST <- c(g,b1,b2)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "NBH-r")     {
    NBH.r <- RAZIAD::reg.model(x, y, b0=b0, dist="NBH-r", link)
    g  <- NBH.r$Estimated.Parameters[1,]                            # Estimated gamma
    b1 <- NBH.r$Estimated.Parameters[2,]                            # Estimated beta
    b2 <- NBH.r$Estimated.Parameters[3,]                            # Estimated beta
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size= y[2], prob=y[3]);};
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*(digamma(x[1]+x[2])+log(x[3])-digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*((x[2]/(x[3]))-(x[1]/(1-x[3])));};
    p01 = function(x) {x[2]^x[1];};
    p0d.r = function(x) {log(x[2])*((x[2])^x[1])}
    p0d.p = function(x) {x[1]*((x[2])^(x[1]-1))}
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    {
    NBH.r <- RAZIAD::reg.model(x, y, b0=b0, dist="NBH-r", link)
    g  <- NBH.r$Estimated.Parameters[1,]                            # Estimated gamma
    b1 <- NBH.r$Estimated.Parameters[2,]                            # Estimated beta
    b2 <- NBH.r$Estimated.Parameters[3,]                            # Estimated beta
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    ptheta1 = function(y){stats::dnbinom(y[1], size= y[2], prob=y[3]);};
    pthetad.r = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*(digamma(x[1]+x[2])+log(x[3])-digamma(x[2]));};
    pthetad.p = function(x){stats::dnbinom(x[1], size = x[2], prob = x[3])*((x[2]/(x[3]))-(x[1]/(1-x[3])));};
    p01 = function(x) {x[2]^x[1];};
    p0d.r = function(x) {log(x[2])*((x[2])^x[1])}
    p0d.p = function(x) {x[1]*((x[2])^(x[1]-1))}
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N));
    xg = x%*%g;                     # X^T*gamma
    xb1 = x1%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    phii = ginv(xg);                # g-inverse
    theta1i = h1inv(xb1)            # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    thetai0 = cbind(0, theta1i, theta2i)
    p00 = apply(thetai0, 1, p01)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x1
    Ci2 = as.vector(h2invd(xb2)) * x

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (log(theta2i))
    A12 = matrix(0, nrow = dim(x)[2], ncol = dim(x)[2])
    A13 = A12
    #In this part we calculate B11
    ##Fisher Information Usingthe other definition and applying trigamma free Method
    #### self-defined function to calculate E Psi_1(nu + Y) given Y ~ NB(r, p)
    M = m
    fEPsi2BNB <- function(nu, r, p, M) {
      trigamma(nu) - sum(stats::pnbinom(0:M, size = r, prob = p, lower.tail = F)/(nu + (0:M))^2);
    }
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = E(y_i + r_i)
    for(i in 1:n){
      nu1 = theta1i[i]
      E11[i] = fEPsi2BNB(nu = nu1, r = theta1i[i], p = theta2i[i], M = M)
    }
    #B11: r*r
    b1b1.0 =  -(E11 + trigamma(theta1i)) * ((1-phii)/(1-p00))
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (p00/(1-p00)) * (1-phii) * p00 * (log(theta1i))^2
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #B12: r*p
    b1b2.0 = -(1/theta2i) * ((1-phii)/(1-p00))
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M1 = (p00/(1-p00)) * (1-phii) * p00 * ((log(theta1i)) * (theta1i/theta2i))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M1))
    B12 = b1b2.1 - M2
    #B22: p*p
    b1b4.0 = (theta1i/((theta2i)^2 * (1-theta2i))) * ((1-phii)/(1-p00))
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (p00/(1-p00)) * (1-phii) * p00 * ((theta1i/theta2i)^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2

    # Fisher Information Matrix
    F1 = cbind(  A11,    A12,    A13)
    F2 = cbind(t(A12),   B11,    B12)
    F3 = cbind(t(A13), t(B12),   B22)
    Fisher = round(rbind(F1, F2, F3), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }

    EST <- c(g,b1,b2)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    phii = ginv(xg);                # g-inverse
    theta1i = h1inv(xb1)            # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    thetai0 = cbind(0, theta1i, theta2i)
    p00 = apply(thetai0, 1, p01)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (log(theta2i))
    A12 = matrix(0, nrow = dim(x)[2], ncol = dim(x)[2])
    A13 = A12
    #In this part we calculate B11
    ##Fisher Information Usingthe other definition and applying trigamma free Method
    #### self-defined function to calculate E Psi_1(nu + Y) given Y ~ NB(r, p)
    M = m
    fEPsi2BNB <- function(nu, r, p, M) {
      trigamma(nu) - sum(stats::pnbinom(0:M, size = r, prob = p, lower.tail = F)/(nu + (0:M))^2);
    }
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = E(y_i + r_i)
    for(i in 1:n){
      nu1 = theta1i[i]
      E11[i] = fEPsi2BNB(nu = nu1, r = theta1i[i], p = theta2i[i], M = M)
    }
    #B11: r*r
    b1b1.0 =  -(E11 + trigamma(theta1i)) * ((1-phii)/(1-p00))
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (p00/(1-p00)) * (1-phii) * p00 * (log(theta1i))^2
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #B12: r*p
    b1b2.0 = -(1/theta2i) * ((1-phii)/(1-p00))
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M1 = (p00/(1-p00)) * (1-phii) * p00 * ((log(theta1i)) * (theta1i/theta2i))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M1))
    B12 = b1b2.1 - M2
    #B22: p*p
    b1b4.0 = (theta1i/((theta2i)^2 * (1-theta2i))) * ((1-phii)/(1-p00))
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (p00/(1-p00)) * (1-phii) * p00 * ((theta1i/theta2i)^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2

    # Fisher Information Matrix
    F1 = cbind(  A11,    A12,    A13)
    F2 = cbind(t(A12),   B11,    B12)
    F3 = cbind(t(A13), t(B12),   B22)
    Fisher = round(rbind(F1, F2, F3), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }

    EST <- c(g,b1,b2)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "ZIBNB")     {
    ZIBNB <- RAZIAD::reg.model(x, y, b0=b0, dist="ZIBNB", link)
    g <- ZIBNB$Estimated.Parameters[1,]                           # Estimated gamma
    b1 <- ZIBNB$Estimated.Parameters[2,]                          # Estimated beta1
    b2 <- ZIBNB$Estimated.Parameters[3,]                          # Estimated beta2
    b3 <- ZIBNB$Estimated.Parameters[4,]                          # Estimated beta3
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h3inv = function(theta) {exp(theta);};                        # h3-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};                       # h3-inverse derivative
    ptheta2 = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]);};        #c(y, n, a,b)
    pthetad.r = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]) * (digamma(x[2]+x[1]) + digamma(x[2]+x[3]) - digamma(x[2]) - digamma(x[1]+x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]) * (digamma(x[2]+x[3]) + digamma(x[3]+x[4]) - digamma(x[1]+x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]) * (digamma(x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[1]+x[3]+x[2]+x[4]) - digamma(x[4]));};
    p02 = function(x){extraDistr::dbnbinom(0, size = x[2], alpha = x[3], beta = x[4]);};               #c(0, n, a,b)
    p0d.r = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) + digamma(x[3]+x[4]) - digamma(x[2]+x[3]+x[4]) + digamma(x[3]));};
    p0d.b = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    xb3 = x%*%b3;                   # X^T*beta3
    phii = ginv(xg);                # g-inverse
    theta1i = h1inv(xb1)            # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    theta3i = h3inv(xb3)            # theta_3 i
    thetai0 = cbind(0, theta1i, theta2i, theta3i)
    p00 = apply(thetai0, 1, p02)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x
    Ci3 = as.vector(h3invd(xb3)) * x

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (digamma(theta1i+theta2i) - digamma(theta1i+theta2i+theta3i))
    A12 = t(A) %*% (Ci1*as.vector(gb10))
    gb20 = (1-phii) * (p00/P) * (digamma(theta1i+theta2i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta2i))
    A13 = t(A) %*% (Ci2 * as.vector(gb20))
    gb30 = (1-phii) * (p00/P) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))
    A14 = t(A) %*% (Ci3 * as.vector(gb30))

    #In this part we calculate B11, B12, B13, B22, B23, B33
    ##Fisher Information Usingthe other definition and applying trigamma free Method
    #### self-defined function to calculate E Psi_1(nu + Y) given Y ~ BNB(r, alpha, beta)
    M = m
    fEPsi3BNB <- function(nu, r, alpha, beta, M) {
      trigamma(nu) - sum(extraDistr::pbnbinom(0:M, size = r, alpha = alpha, beta = beta, lower.tail = F)/(nu + (0:M))^2);
    }
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = E(y_i + r_i + alpha_i + beta_i)
    E22 = matrix(NA, nrow = n, ncol = 1) # E22 = E(y_i + beta_i)
    E33 = matrix(NA, nrow = n, ncol=1)   # E33 = E(y_i + r_i)
    for(i in 1:n){
      nu1 = theta1i[i] + theta2i[i] + theta3i[i]
      E11[i] = fEPsi3BNB(nu = nu1, r = theta1i[i], alpha = theta2i[i], beta = theta3i[i], M = M)
      nu2 = theta3i[i]
      E22[i] = fEPsi3BNB(nu = nu2, r = theta1i[i], alpha = theta2i[i], beta = theta3i[i], M = M)
      nu3 = theta1i[i]
      E33[i] = fEPsi3BNB(nu = nu3, r = theta1i[i], alpha = theta2i[i], beta = theta3i[i], M = M)
    }

    #B11: r*r
    b1b1.0 =  -(E33 + trigamma(theta1i + theta2i) - trigamma(theta1i) - E11) * (1-phii)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta2i) - digamma(theta1i+theta2i+theta3i))^2)
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #B12: r*a
    b1b2.0 = -(trigamma(theta1i + theta2i) - E11) * (1-phii)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M1 = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta2i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta1i+theta2i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta2i)))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M1))
    B12 = b1b2.1 - M2
    #B13: r*b
    b1b3.0 = E11 * (1-phii)
    b1b3.1 = t(Ci1) %*% (Ci3 * as.vector(b1b3.0))
    N = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta2i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    N2 = t(Ci1) %*% (Ci3 * as.vector(N))
    B13 = b1b3.1 - N2
    #B22: a*a
    b1b4.0 = -(trigamma(theta1i + theta2i) - E11  + trigamma(theta2i + theta3i) - trigamma(theta2i)) * (1-phii)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta2i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta2i))^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2
    #B23: a*b
    b1b5.0 = -(-E11 + trigamma(theta2i + theta3i)) * (1-phii)
    b1b5.1 = t(Ci2) %*% (Ci3 * as.vector(b1b5.0))
    O = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta2i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta2i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    O2 = t(Ci2) %*% (Ci3 * as.vector(O))
    B23 = b1b5.1 - O2
    #B33: b*b
    b2b6.0 = -(E22 - E11 + trigamma(theta2i + theta3i) - trigamma(theta3i)) * (1-phii)
    b2b6.1 = t(Ci3) %*% (Ci3 * as.vector(b2b6.0))
    Q = (1/P) * phii * (1-phii) * p00 * ((digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    Q2 = t(Ci3) %*% (Ci3 * as.vector(Q))
    B44 = b2b6.1 - Q2

    # Fisher Information Matrix
    F1 = cbind(  A11,    A12,    A13,   A14)
    F2 = cbind(t(A12),   B11,    B12,   B13)
    F3 = cbind(t(A13), t(B12),   B22,   B23)
    F4 = cbind(t(A14), t(B13),  t(B23), B44)
    Fisher = round(rbind(F1, F2, F3, F4), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }

    EST <- c(g,b1,b2,b3)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "BNBH")      {
    BNBH <- RAZIAD::reg.model(x, y, b0=b0, dist="BNBH", link)
    g <- BNBH$Estimated.Parameters[1,]                            # Estimated gamma
    b1 <- BNBH$Estimated.Parameters[2,]                           # Estimated beta1
    b2 <- BNBH$Estimated.Parameters[3,]                           # Estimated beta2
    b3 <- BNBH$Estimated.Parameters[4,]                           # Estimated beta3
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h3inv = function(theta) {exp(theta);};                        # h3-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};                       # h3-inverse derivative
    ptheta2 = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]);};        #c(y, n, a,b)
    pthetad.r = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]) * (digamma(x[2]+x[1]) + digamma(x[2]+x[3]) - digamma(x[2]) - digamma(x[1]+x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]) * (digamma(x[2]+x[3]) + digamma(x[3]+x[4]) - digamma(x[1]+x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbnbinom(x[1], size = x[2], alpha = x[3], beta = x[4]) * (digamma(x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[1]+x[3]+x[2]+x[4]) - digamma(x[4]));};
    p02 = function(x){extraDistr::dbnbinom(0, size = x[2], alpha = x[3], beta = x[4]);};               #c(0, n, a,b)
    p0d.r = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[2]+x[3]) + digamma(x[3]+x[4]) - digamma(x[2]+x[3]+x[4]) + digamma(x[3]));};
    p0d.b = function(x){extraDistr::dbnbinom(0, size=x[2], alpha=x[3], beta=x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    xb3 = x%*%b3;                   # X^T*beta3
    phii = ginv(xg);                # g-inverse
    theta1i = h1inv(xb1)            # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    theta3i = h3inv(xb3)            # theta_3 i
    thetai0 = cbind(0, theta1i, theta2i, theta3i)
    p00 = apply(thetai0, 1, p02)

    P1 = phii * (1 - phii)
    P2 = (1 - phii)/(1 - p00)
    A = as.vector(ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x
    Ci3 = as.vector(h3invd(xb3)) * x

    A11 = t(A) %*% (A * as.vector(1/P1))
    A12 = matrix(0, nrow = dim(x)[2], ncol = dim(x)[2])
    A13 = A12
    A14 = A12

    #In this part we calculate B11, B12, B13, B22, B23, B33
    ##Fisher Information Usingthe other definition and applying trigamma free Method
    #### self-defined function to calculate E Psi_1(nu + Y) given Y ~ BNB(r, alpha, beta)
    M = m
    fEPsi4BNB <- function(nu, r, alpha, beta, M) {
      trigamma(nu) - sum(extraDistr::pbnbinom(0:M, size = r, alpha = alpha, beta = beta, lower.tail = F)/(nu + (0:M))^2);
    }
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = E(y_i + r_i + alpha_i + beta_i)
    E22 = matrix(NA, nrow = n, ncol = 1) # E22 = E(y_i + beta_i)
    E33 = matrix(NA, nrow = n, ncol=1)   # E33 = E(y_i + r_i)
    for(i in 1:n){
      nu1 = theta1i[i] + theta2i[i] + theta3i[i]
      E11[i] = fEPsi4BNB(nu = nu1, r = theta1i[i], alpha = theta2i[i], beta = theta3i[i], M = M)
      nu2 = theta3i[i]
      E22[i] = fEPsi4BNB(nu = nu2, r = theta1i[i], alpha = theta2i[i], beta = theta3i[i], M = M)
      nu3 = theta1i[i]
      E33[i] = fEPsi4BNB(nu = nu3, r = theta1i[i], alpha = theta2i[i], beta = theta3i[i], M = M)
    }

    #B11: r*r
    b1b1.0 =  -(E33 + trigamma(theta1i + theta2i) - trigamma(theta1i) - E11) * (P2)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta2i) - digamma(theta1i+theta2i+theta3i))^2)
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #B12: r*a
    b1b2.0 = -(trigamma(theta1i + theta2i) - E11) * (P2)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M1 = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta2i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta1i+theta2i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta2i)))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M1))
    B12 = b1b2.1 - M2
    #B13: r*b
    b1b3.0 = E11 * (P2)
    b1b3.1 = t(Ci1) %*% (Ci3 * as.vector(b1b3.0))
    N = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta2i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    N2 = t(Ci1) %*% (Ci3 * as.vector(N))
    B13 = b1b3.1 - N2
    #B22: a*a
    b1b4.0 = -(trigamma(theta1i + theta2i) - E11  + trigamma(theta2i + theta3i) - trigamma(theta2i)) * (P2)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    R = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta2i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta2i))^2)
    R2 = t(Ci2) %*% (Ci2 * as.vector(R))
    B22 = b1b4.1 - R2
    #B23: a*b
    b1b5.0 = -(-E11 + trigamma(theta2i + theta3i)) * (P2)
    b1b5.1 = t(Ci2) %*% (Ci3 * as.vector(b1b5.0))
    O = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta2i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta2i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    O2 = t(Ci2) %*% (Ci3 * as.vector(O))
    B23 = b1b5.1 - O2
    #B33: b*b
    b2b6.0 = -(E22 - E11 + trigamma(theta2i + theta3i) - trigamma(theta3i)) * (P2)
    b2b6.1 = t(Ci3) %*% (Ci3 * as.vector(b2b6.0))
    Q = (p00/(1-p00)) * (P2) * ((digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    Q2 = t(Ci3) %*% (Ci3 * as.vector(Q))
    B44 = b2b6.1 - Q2

    # Fisher Information Matrix
    F1 = cbind(  A11,    A12,    A13,   A14)
    F2 = cbind(matrix(A12, nrow = dim(x)[2]),   B11,    B12,   B13)
    F3 = cbind(matrix(A13, nrow = dim(x)[2]), t(B12),   B22,   B23)
    F4 = cbind(matrix(A14, nrow = dim(x)[2]), t(B13),  t(B23), B44)
    Fisher = round(rbind(F1, F2, F3, F4), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }

    EST <- c(g,b1,b2,b3)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "ZIBB")      {
    ZIBB <- RAZIAD::reg.model(x, y, b0=b0, dist="ZIBB", link)
    g <- ZIBB$Estimated.Parameters[1,]                            # Estimated gamma
    b1 <- ZIBB$Estimated.Parameters[2,]                           # Estimated beta1
    b2 <- ZIBB$Estimated.Parameters[3,]                           # Estimated beta2
    b3 <- ZIBB$Estimated.Parameters[4,]                           # Estimated beta3
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h3inv = function(theta) {exp(theta);};                        # h3-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};                       # h3-inverse derivative
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};        #c(y, n, a,b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];                 # number of predictors
    N = dim(x)[1];
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    xb3 = x%*%b3;                   # X^T*beta3
    phii = ginv(xg);                # g-inverse
    theta1i = ceiling(h1inv(xb1))   # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    theta3i = h3inv(xb3)            # theta_3 i
    thetai0 = cbind(0, theta1i, theta2i, theta3i)
    p00 = apply(thetai0, 1, p03)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x
    Ci3 = as.vector(h3invd(xb3)) * x

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))
    A12 = t(A) %*% (Ci1*as.vector(gb10))
    gb20 = (1-phii) * (p00/P) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))
    A13 = t(A) %*% (Ci2 * as.vector(gb20))
    gb30 = (1-phii) * (p00/P) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))
    A14 = t(A) %*% (Ci3 * as.vector(gb30))

    fEPsi1BB1 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(nu - 0:n))}
    fEPsi1BB2 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(0:n + nu))}
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = EPsi_1(n_i + 1 - y_i)
    E22 = matrix(NA, nrow = n, ncol = 1) # E22 = EPsi_1(n_i + b_i - y_i)
    E33 = matrix(NA, nrow = n, ncol=1)   # E33 = EPsi_1(y_i + alpha_i)
    for(i in 1:n){
      E11[i] = fEPsi1BB1(nu = theta1i[i]+1, n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])          #EPsi_1(n_i + 1 - y_i)
      E22[i] = fEPsi1BB1(nu = theta1i[i]+theta3i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i]) #E22 = EPsi_1(n_i + b_i - y_i)
      E33[i] = fEPsi1BB2(nu = theta2i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])            #E33 = EPsi_1(y_i + alpha_i)
    }

    #n*n
    b1b1.0 = (-trigamma(theta1i + 1) + E11 - E22 + trigamma(theta1i + theta2i + theta2i)) * (1-phii)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #n*a
    b1b2.0 = trigamma(theta1i + theta2i + theta3i) * (1-phii)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M))
    B12 = b1b2.1 - M2
    #n*b
    b1b3.0 = (-E22 + trigamma(theta1i + theta2i + theta3i)) * (1-phii)
    b1b3.1 = t(Ci1) %*% (Ci3 * as.vector(b1b3.0))
    N = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    N2 = t(Ci1) %*% (Ci3 * as.vector(N))
    B13 = b1b3.1 - N2
    #a*a
    b1b4.0 = (-E33 + trigamma(theta1i + theta2i + theta3i) -trigamma(theta2i + theta3i) + trigamma(theta2i)) * (1-phii)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (1/P) * phii * (1-phii) * p00 * ((digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2
    #a*b
    b1b5.0 = (trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i)) * (1-phii)
    b1b5.1 = t(Ci2) %*% (Ci3 * as.vector(b1b5.0))
    O = (1/P) * phii * (1-phii) * p00 * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    O2 = t(Ci2) %*% (Ci3 * as.vector(O))
    B23 = b1b5.1 - O2
    #b*b
    b2b6.0 = (-E22 + trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i) + trigamma(theta3i)) * (1-phii)
    b2b6.1 = t(Ci3) %*% (Ci3 * as.vector(b2b6.0))
    Q = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))^2)
    Q2 = t(Ci3) %*% (Ci3 * as.vector(Q))
    B44 = b2b6.1 - Q2

    F1 = cbind(  A11,    A12,    A13,   A14)
    F2 = cbind(t(A12),   B11,    B12,   B13)
    F3 = cbind(t(A13), t(B12),   B22,   B23)
    F4 = cbind(t(A14), t(B13),  t(B23), B44)

    # Check if Fisher is positive definite
    Fisher = round(rbind(F1, F2, F3, F4),2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }
    EST <- c(g,b1,b2,b3)
    df = length(EST)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt <- round(cbind(Conf, EST, Sig, S_E, z_score),4)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))}
  if (dist == "BBH")       {
    BBH <- RAZIAD::reg.model(x, y, b0=b0, dist="BBH", link)
    g <- BBH$Estimated.Parameters[1,]                             # Estimated gamma
    b1 <- BBH$Estimated.Parameters[2,]                            # Estimated beta1
    b2 <- BBH$Estimated.Parameters[3,]                            # Estimated beta2
    b3 <- BBH$Estimated.Parameters[4,]                            # Estimated beta3
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h3inv = function(theta) {exp(theta);};                        # h3-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};                       # h3-inverse derivative
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};        #c(y, n, a,b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    #p0=function(x) {gamma(x[2]+x[4])*gamma(x[3]+x[4])/(gamma(x[2]+x[3]+x[4])*gamma(x[4]));}; #x[1]=n
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                   # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    xb3 = x%*%b3;                   # X^T*beta3
    phii = ginv(xg);                # g-inverse
    theta1i = ceiling(h1inv(xb1))   # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    theta3i = h3inv(xb3)            # theta_3 i
    thetai0 = cbind(0, theta1i, theta2i, theta3i)
    p00 = apply(thetai0, 1, p03)

    P1 = phii * (1 - phii)
    P2 = (1 - phii)/(1 - p00)
    A = as.vector(ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x
    Ci3 = as.vector(h3invd(xb3)) * x

    A11 = t(A) %*% (A * as.vector(1/P1))
    A12 = matrix(0, nrow = dim(x)[2], ncol = dim(x)[2])
    A13 = A12
    A14 = A12

    fEPsi1BB1 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(nu - 0:n))}
    fEPsi1BB2 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(0:n + nu))}
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = EPsi_1(n_i + 1 - y_i)
    E22 = matrix(NA, nrow = n, ncol = 1) # E22 = EPsi_1(n_i + b_i - y_i)
    E33 = matrix(NA, nrow = n, ncol=1)   # E33 = EPsi_1(y_i + alpha_i)
    for(i in 1:n){
      E11[i] = fEPsi1BB1(nu = theta1i[i]+1, n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])          #EPsi_1(n_i + 1 - y_i)
      E22[i] = fEPsi1BB1(nu = theta1i[i]+theta3i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i]) #E22 = EPsi_1(n_i + b_i - y_i)
      E33[i] = fEPsi1BB2(nu = theta2i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])            #E33 = EPsi_1(y_i + alpha_i)
    }

    #n*n
    b1b1.0 = (-trigamma(theta1i + 1) + E11 - E22 + trigamma(theta1i + theta2i + theta2i)) * (P2)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #n*a
    b1b2.0 = trigamma(theta1i + theta2i + theta3i) * (P2)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M))
    B12 = b1b2.1 - M2
    #n*b
    b1b3.0 = (-E22 + trigamma(theta1i + theta2i + theta3i)) * (P2)
    b1b3.1 = t(Ci1) %*% (Ci3 * as.vector(b1b3.0))
    N = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    N2 = t(Ci1) %*% (Ci3 * as.vector(N))
    B13 = b1b3.1 - N2
    #a*a
    b1b4.0 = (-E33 + trigamma(theta1i + theta2i + theta3i) -trigamma(theta2i + theta3i) + trigamma(theta2i)) * (P2)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    R = (p00/(1-p00)) * (P2) * ((digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    R2 = t(Ci2) %*% (Ci2 * as.vector(R))
    B22 = b1b4.1 - R2
    #a*b
    b1b5.0 = (trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i)) * (P2)
    b1b5.1 = t(Ci2) %*% (Ci3 * as.vector(b1b5.0))
    O = (p00/(1-p00)) * (P2) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    O2 = t(Ci2) %*% (Ci3 * as.vector(O))
    B23 = b1b5.1 - O2
    #b*b
    b2b6.0 = (-E22 + trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i) + trigamma(theta3i)) * (P2)
    b2b6.1 = t(Ci3) %*% (Ci3 * as.vector(b2b6.0))
    Q = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))^2)
    Q2 = t(Ci3) %*% (Ci3 * as.vector(Q))
    B44 = b2b6.1 - Q2

    F1 = cbind(  A11,    A12,    A13,   A14)
    F2 = cbind(t(A12),   B11,    B12,   B13)
    F3 = cbind(t(A13), t(B12),   B22,   B23)
    F4 = cbind(t(A14), t(B13),  t(B23), B44)

    # Check if Fisher is positive definite
    Fisher=round(rbind(F1, F2, F3, F4), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }
    EST <- c(g,b1,b2,b3)
    df = length(EST)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*4), ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E
    ConfInt <- round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))}
  if (dist == "ZIBB-n")    {
    ZIBB.n <- RAZIAD::reg.model(x, y, b0=b0, dist="ZIBB-n", link)
    g <- ZIBB.n$Estimated.Gamma                                   # Estimated gamma
    b1 <- ZIBB.n$Estimated.Beta1                                  # Estimated beta1
    b2 <- ZIBB.n$Estimated.Beta2                                  # Estimated beta2
    b3 <- ZIBB.n$Estimated.Beta3                                  # Estimated beta3
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h3inv = function(theta) {exp(theta);};                        # h3-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};                       # h3-inverse derivative
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};        #c(y, n, a,b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    #p0=function(x) {gamma(x[2]+x[4])*gamma(x[3]+x[4])/(gamma(x[2]+x[3]+x[4])*gamma(x[4]));}; #x[1]=n
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N));
    xg = x%*%g;                     # X^T*gamma
    xb1 = x1%*%b1;                  # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    xb3 = x%*%b3;                   # X^T*beta3
    phii = ginv(xg);                # g-inverse
    theta1i = ceiling(h1inv(xb1))   # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    theta3i = h3inv(xb3)            # theta_3 i
    thetai0 = cbind(0, theta1i, theta2i, theta3i)
    p00 = apply(thetai0, 1, p03)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x1
    Ci2 = as.vector(h2invd(xb2)) * x
    Ci3 = as.vector(h3invd(xb3)) * x

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))
    A12 = t(A) %*% (Ci1*as.vector(gb10))
    gb20 = (1-phii) * (p00/P) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))
    A13 = t(A) %*% (Ci2 * as.vector(gb20))
    gb30 = (1-phii) * (p00/P) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))
    A14 = t(A) %*% (Ci3 * as.vector(gb30))

    fEPsi1BB1 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(nu - 0:n))}
    fEPsi1BB2 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(0:n + nu))}
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = EPsi_1(n_i + 1 - y_i)
    E22 = matrix(NA, nrow = n, ncol = 1) # E22 = EPsi_1(n_i + b_i - y_i)
    E33 = matrix(NA, nrow = n, ncol=1)   # E33 = EPsi_1(y_i + alpha_i)
    for(i in 1:n){
      E11[i] = fEPsi1BB1(nu = theta1i[i]+1, n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])          #EPsi_1(n_i + 1 - y_i)
      E22[i] = fEPsi1BB1(nu = theta1i[i]+theta3i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i]) #E22 = EPsi_1(n_i + b_i - y_i)
      E33[i] = fEPsi1BB2(nu = theta2i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])            #E33 = EPsi_1(y_i + alpha_i)
    }

    #n*n
    b1b1.0 = (-trigamma(theta1i + 1) + E11 - E22 + trigamma(theta1i + theta2i + theta2i)) * (1-phii)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #n*a
    b1b2.0 = trigamma(theta1i + theta2i + theta3i) * (1-phii)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M))
    B12 = b1b2.1 - M2
    #n*b
    b1b3.0 = (-E22 + trigamma(theta1i + theta2i + theta3i)) * (1-phii)
    b1b3.1 = t(Ci1) %*% (Ci3 * as.vector(b1b3.0))
    N = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    N2 = t(Ci1) %*% (Ci3 * as.vector(N))
    B13 = b1b3.1 - N2
    #a*a
    b1b4.0 = (-E33 + trigamma(theta1i + theta2i + theta3i) -trigamma(theta2i + theta3i) + trigamma(theta2i)) * (1-phii)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (1/P) * phii * (1-phii) * p00 * ((digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2
    #a*b
    b1b5.0 = (trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i)) * (1-phii)
    b1b5.1 = t(Ci2) %*% (Ci3 * as.vector(b1b5.0))
    O = (1/P) * phii * (1-phii) * p00 * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    O2 = t(Ci2) %*% (Ci3 * as.vector(O))
    B23 = b1b5.1 - O2
    #b*b
    b2b6.0 = (-E22 + trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i) + trigamma(theta3i)) * (1-phii)
    b2b6.1 = t(Ci3) %*% (Ci3 * as.vector(b2b6.0))
    Q = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))^2)
    Q2 = t(Ci3) %*% (Ci3 * as.vector(Q))
    B44 = b2b6.1 - Q2

    F1 = cbind(  A11,    A12,    A13,   A14)
    F2 = cbind(t(A12),   B11,    B12,   B13)
    F3 = cbind(t(A13), t(B12),   B22,   B23)
    F4 = cbind(t(A14), t(B13),  t(B23), B44)

    Fisher=round(rbind(F1, F2, F3, F4), 2)
    # is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)
    #
    # if (is_positive_definite) {
    #   # If Fisher is already positive definite, take its inverse
    #   Fisher_inv <- solve(Fisher)
    # } else {
    #   # If Fisher is not positive definite, make it positive definite
    #   Fisher_pos_def <- Matrix::nearPD(Fisher)
    #   Fisher_pd <- Fisher_pos_def$mat
    #   Fisher_inv <- solve(Fisher_pd)
    # }
    Fisher_inv <- solve(Fisher)
    EST<-c(g,b1,b2,b3)
    df = length(EST)
    A1 = EST + 1.96 * sqrt(diag(abs(Fisher_inv)))
    A2 = EST - 1.96 * sqrt(diag(abs(Fisher_inv)))
    Conf = matrix(c(A2,A1), nrow = (np*3)+1, ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(abs(Fisher_inv)))
    z_score <- EST / S_E  # Calculate z-score
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "BBH-n")     {
    BBH.n <- RAZIAD::reg.model(x, y, b0=b0, dist="BBH-n", link)
    g <- BBH.n$Estimated.Gamma                                    # Estimated gamma
    b1 <- BBH.n$Estimated.Beta1                                   # Estimated beta1
    b2 <- BBH.n$Estimated.Beta2                                   # Estimated beta2
    b3 <- BBH.n$Estimated.Beta3                                   # Estimated beta3
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h3inv = function(theta) {exp(theta);};                        # h3-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};                       # h3-inverse derivative
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};        #c(y, n, a,b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    #p0=function(x) {gamma(x[2]+x[4])*gamma(x[3]+x[4])/(gamma(x[2]+x[3]+x[4])*gamma(x[4]));}; #x[1]=n
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N));
    xg = x%*%g;                     # X^T*gamma
    xb1 = x1%*%b1;                  # X^T*beta1
    xb2 = x%*%b2;                   # X^T*beta2
    xb3 = x%*%b3;                   # X^T*beta3
    phii = ginv(xg);                # g-inverse
    theta1i = ceiling(h1inv(xb1))   # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    theta3i = h3inv(xb3)            # theta_3 i
    thetai0 = cbind(0, theta1i, theta2i, theta3i)
    p00 = apply(thetai0, 1, p03)

    P1 = phii * (1 - phii)
    P2 = (1 - phii)/(1 - p00)
    A = as.vector(ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x1
    Ci2 = as.vector(h2invd(xb2)) * x
    Ci3 = as.vector(h3invd(xb3)) * x

    A11 = t(A) %*% (A * as.vector(1/P1))
    A12 = matrix(0, nrow = dim(x)[2], ncol = dim(x)[2])
    A13 = A12
    A14 = A12

    fEPsi1BB1 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(nu - 0:n))}
    fEPsi1BB2 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(0:n + nu))}
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = EPsi_1(n_i + 1 - y_i)
    E22 = matrix(NA, nrow = n, ncol = 1) # E22 = EPsi_1(n_i + b_i - y_i)
    E33 = matrix(NA, nrow = n, ncol=1)   # E33 = EPsi_1(y_i + alpha_i)
    for(i in 1:n){
      E11[i] = fEPsi1BB1(nu = theta1i[i]+1, n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])          #EPsi_1(n_i + 1 - y_i)
      E22[i] = fEPsi1BB1(nu = theta1i[i]+theta3i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i]) #E22 = EPsi_1(n_i + b_i - y_i)
      E33[i] = fEPsi1BB2(nu = theta2i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])            #E33 = EPsi_1(y_i + alpha_i)
    }

    #n*n
    b1b1.0 = (-trigamma(theta1i + 1) + E11 - E22 + trigamma(theta1i + theta2i + theta2i)) * (P2)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #n*a
    b1b2.0 = trigamma(theta1i + theta2i + theta3i) * (P2)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M))
    B12 = b1b2.1 - M2
    #n*b
    b1b3.0 = (-E22 + trigamma(theta1i + theta2i + theta3i)) *  (P2)
    b1b3.1 = t(Ci1) %*% (Ci3 * as.vector(b1b3.0))
    N = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    N2 = t(Ci1) %*% (Ci3 * as.vector(N))
    B13 = b1b3.1 - N2
    #a*a
    b1b4.0 = (-E33 + trigamma(theta1i + theta2i + theta3i) -trigamma(theta2i + theta3i) + trigamma(theta2i)) *  (P2)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    R = (p00/(1-p00)) * (P2) * ((digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    R2 = t(Ci2) %*% (Ci2 * as.vector(R))
    B22 = b1b4.1 - R2
    #a*b
    b1b5.0 = (trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i)) *  (P2)
    b1b5.1 = t(Ci2) %*% (Ci3 * as.vector(b1b5.0))
    O = (p00/(1-p00)) * (P2) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    O2 = t(Ci2) %*% (Ci3 * as.vector(O))
    B23 = b1b5.1 - O2
    #b*b
    b2b6.0 = (-E22 + trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i) + trigamma(theta3i)) *  (P2)
    b2b6.1 = t(Ci3) %*% (Ci3 * as.vector(b2b6.0))
    Q = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))^2)
    Q2 = t(Ci3) %*% (Ci3 * as.vector(Q))
    B44 = b2b6.1 - Q2

    F1 = cbind(  A11,    A12,    A13,   A14)
    F2 = cbind(t(A12),   B11,    B12,   B13)
    F3 = cbind(t(A13), t(B12),   B22,   B23)
    F4 = cbind(t(A14), t(B13),  t(B23), B44)

    Fisher=round(rbind(F1, F2, F3, F4), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }
    EST<-c(g,b1,b2,b3)
    df = length(EST)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*3)+1, ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv))
    z_score <- EST / S_E  # Calculate z-score
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "ZIBB-ab")   {
    ZIBB.n <- RAZIAD::reg.model(x, y, b0=b0, dist="ZIBB-ab", link)
    g <- ZIBB.n$Estimated.Gamma                                   # Estimated gamma
    b1 <- ZIBB.n$Estimated.Beta1                                  # Estimated beta1
    b2 <- ZIBB.n$Estimated.Beta2                                  # Estimated beta2
    b3 <- ZIBB.n$Estimated.Beta3                                  # Estimated beta3
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h3inv = function(theta) {exp(theta);};                        # h3-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};                       # h3-inverse derivative
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};        #c(y, n, a,b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    #p0=function(x) {gamma(x[2]+x[4])*gamma(x[3]+x[4])/(gamma(x[2]+x[3]+x[4])*gamma(x[4]));}; #x[1]=n
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N));
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                  # X^T*beta1
    xb2 = x1%*%b2;                   # X^T*beta2
    xb3 = x1%*%b3;                   # X^T*beta3
    phii = ginv(xg);                # g-inverse
    theta1i = ceiling(h1inv(xb1))   # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    theta3i = h3inv(xb3)            # theta_3 i
    thetai0 = cbind(0, theta1i, theta2i, theta3i)
    p00 = apply(thetai0, 1, p03)

    P = phii + (1 - phii) * p00
    A = as.vector((1/(1 - phii)) * ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x1
    Ci3 = as.vector(h3invd(xb3)) * x1

    A11 = t(A) %*% (A * as.vector((1/P) - 1))
    gb10 = (1-phii) * (p00/P) * (digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))
    A12 = t(A) %*% (Ci1*as.vector(gb10))
    gb20 = (1-phii) * (p00/P) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))
    A13 = t(A) %*% (Ci2 * as.vector(gb20))
    gb30 = (1-phii) * (p00/P) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))
    A14 = t(A) %*% (Ci3 * as.vector(gb30))

    fEPsi1BB1 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(nu - 0:n))}
    fEPsi1BB2 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(0:n + nu))}
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = EPsi_1(n_i + 1 - y_i)
    E22 = matrix(NA, nrow = n, ncol = 1) # E22 = EPsi_1(n_i + b_i - y_i)
    E33 = matrix(NA, nrow = n, ncol=1)   # E33 = EPsi_1(y_i + alpha_i)
    for(i in 1:n){
      E11[i] = fEPsi1BB1(nu = theta1i[i]+1, n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])          #EPsi_1(n_i + 1 - y_i)
      E22[i] = fEPsi1BB1(nu = theta1i[i]+theta3i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i]) #E22 = EPsi_1(n_i + b_i - y_i)
      E33[i] = fEPsi1BB2(nu = theta2i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i]) #E33 = EPsi_1(y_i + alpha_i)
    }

    #n*n
    b1b1.0 = (-trigamma(theta1i + 1) + E11 - E22 + trigamma(theta1i + theta2i + theta2i)) * (1-phii)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #n*a
    b1b2.0 = trigamma(theta1i + theta2i + theta3i) * (1-phii)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M))
    B12 = b1b2.1 - M2
    #n*b
    b1b3.0 = (-E22 + trigamma(theta1i + theta2i + theta3i)) * (1-phii)
    b1b3.1 = t(Ci1) %*% (Ci3 * as.vector(b1b3.0))
    N = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    N2 = t(Ci1) %*% (Ci3 * as.vector(N))
    B13 = b1b3.1 - N2
    #a*a
    b1b4.0 = (-E33 + trigamma(theta1i + theta2i + theta3i) -trigamma(theta2i + theta3i) + trigamma(theta2i)) * (1-phii)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    P1 = (1/P) * phii * (1-phii) * p00 * ((digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    P2 = t(Ci2) %*% (Ci2 * as.vector(P1))
    B22 = b1b4.1 - P2
    #a*b
    b1b5.0 = (trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i)) * (1-phii)
    b1b5.1 = t(Ci2) %*% (Ci3 * as.vector(b1b5.0))
    O = (1/P) * phii * (1-phii) * p00 * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    O2 = t(Ci2) %*% (Ci3 * as.vector(O))
    B23 = b1b5.1 - O2
    #b*b
    b2b6.0 = (-E22 + trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i) + trigamma(theta3i)) * (1-phii)
    b2b6.1 = t(Ci3) %*% (Ci3 * as.vector(b2b6.0))
    Q = (1/P) * phii * (1-phii) * p00 * ((digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))^2)
    Q2 = t(Ci3) %*% (Ci3 * as.vector(Q))
    B44 = b2b6.1 - Q2

    F1 = cbind(  A11,    A12,    A13,   A14)
    F2 = cbind(t(A12),   B11,    B12,   B13)
    F3 = cbind(t(A13), t(B12),   B22,   B23)
    F4 = cbind(t(A14), t(B13),  t(B23), B44)

    Fisher=round(rbind(F1, F2, F3, F4), 2)
    Fisher=round(rbind(F1, F2, F3, F4), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }
    EST<-c(g,b1,b2,b3)
    df = length(EST)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*2)+2, ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv)/n)
    z_score <- EST / S_E  # Calculate z-score
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
  if (dist == "BBH-ab")    {
    BBH.n <- RAZIAD::reg.model(x, y, b0=b0, dist="BBH-ab", link)
    g <- BBH.n$Estimated.Gamma                                    # Estimated gamma
    b1 <- BBH.n$Estimated.Beta1                                   # Estimated beta1
    b2 <- BBH.n$Estimated.Beta2                                   # Estimated beta2
    b3 <- BBH.n$Estimated.Beta3                                   # Estimated beta3
    h1inv = function(theta) {exp(theta);};                        # h1-inverse
    h2inv = function(theta) {exp(theta);};                        # h2-inverse
    h3inv = function(theta) {exp(theta);};                        # h3-inverse
    h1invd = function(theta) {exp(theta);};                       # h1-inverse derivative
    h2invd = function(theta) {exp(theta);};                       # h2-inverse derivative
    h3invd = function(theta) {exp(theta);};                       # h3-inverse derivative
    ptheta3 = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]);};        #c(y, n, a,b)
    pthetad.n = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+1) + digamma(x[2]-x[1]+x[4]) - digamma(x[2]-x[1]+1) - digamma(x[3]+x[2]+x[4]));};
    pthetad.a = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[1]+x[3]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[3]));};
    pthetad.b = function(x){extraDistr::dbbinom(x[1], size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]-x[1]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    #p0=function(x) {gamma(x[2]+x[4])*gamma(x[3]+x[4])/(gamma(x[2]+x[3]+x[4])*gamma(x[4]));}; #x[1]=n
    p03 = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]);};      # c(0, n, a, b)
    p0d.n = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.a = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]));};
    p0d.b = function(x){extraDistr::dbbinom(0, size = ceiling(x[2]), alpha = x[3], beta = x[4]) * (digamma(x[2]+x[4]) + digamma(x[3]+x[4]) - digamma(x[3]+x[2]+x[4]) - digamma(x[4]));};
    np = dim(x)[2];  # number of predictors
    N = dim(x)[1];
    x1=as.matrix(rep(1, N));
    xg = x%*%g;                     # X^T*gamma
    xb1 = x%*%b1;                   # X^T*beta1
    xb2 = x1%*%b2;                  # X^T*beta2
    xb3 = x1%*%b3;                  # X^T*beta3
    phii = ginv(xg);                # g-inverse
    theta1i = ceiling(h1inv(xb1))   # theta_1 i
    theta2i = h2inv(xb2)            # theta_2 i
    theta3i = h3inv(xb3)            # theta_3 i
    thetai0 = cbind(0, theta1i, theta2i, theta3i)
    p00 = apply(thetai0, 1, p03)

    P1 = phii * (1 - phii)
    P2 = (1 - phii)/(1 - p00)
    A = as.vector(ginvd(xg)) * x
    Ci1 = as.vector(h1invd(xb1)) * x
    Ci2 = as.vector(h2invd(xb2)) * x1
    Ci3 = as.vector(h3invd(xb3)) * x1

    A11 = t(A) %*% (A * as.vector(1/P1))
    A12 = matrix(0, nrow = dim(x)[2], ncol = dim(x)[2])
    A13 = A12
    A14 = A12

    fEPsi1BB1 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(nu - 0:n))}
    fEPsi1BB2 <- function(nu, n, alpha, beta){sum(extraDistr::dbbinom(0:n, n, alpha, beta) * trigamma(0:n + nu))}
    n = dim(x)[1]
    E11 = matrix(NA, nrow = n, ncol = 1) # E11 = EPsi_1(n_i + 1 - y_i)
    E22 = matrix(NA, nrow = n, ncol = 1) # E22 = EPsi_1(n_i + b_i - y_i)
    E33 = matrix(NA, nrow = n, ncol=1)   # E33 = EPsi_1(y_i + alpha_i)
    for(i in 1:n){
      E11[i] = fEPsi1BB1(nu = theta1i[i]+1, n = theta1i[i], alpha = theta2i[i], beta = theta3i[i])          #EPsi_1(n_i + 1 - y_i)
      E22[i] = fEPsi1BB1(nu = theta1i[i]+theta3i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i]) #E22 = EPsi_1(n_i + b_i - y_i)
      E33[i] = fEPsi1BB2(nu = theta2i[i], n = theta1i[i], alpha = theta2i[i], beta = theta3i[i]) #E33 = EPsi_1(y_i + alpha_i)
    }

    #n*n
    b1b1.0 = (-trigamma(theta1i + 1) + E11 - E22 + trigamma(theta1i + theta2i + theta2i)) * (P2)
    b1b1.1 = t(Ci1) %*% (Ci1 * as.vector(b1b1.0))
    H = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    H2 = t(Ci1) %*% (Ci1 * as.vector(H))
    B11 = b1b1.1 - H2
    #n*a
    b1b2.0 = trigamma(theta1i + theta2i + theta3i) * (P2)
    b1b2.1 = t(Ci1) %*% (Ci2 * as.vector(b1b2.0))
    M = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i)))
    M2 = t(Ci1) %*% (Ci2 * as.vector(M))
    B12 = b1b2.1 - M2
    #n*b
    b1b3.0 = (-E22 + trigamma(theta1i + theta2i + theta3i)) * (P2)
    b1b3.1 = t(Ci1) %*% (Ci3 * as.vector(b1b3.0))
    N = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) - digamma(theta1i+theta2i+theta3i)) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    N2 = t(Ci1) %*% (Ci3 * as.vector(N))
    B13 = b1b3.1 - N2
    #a*a
    b1b4.0 = (-E33 + trigamma(theta1i + theta2i + theta3i) -trigamma(theta2i + theta3i) + trigamma(theta2i)) * (P2)
    b1b4.1 = t(Ci2) %*% (Ci2 * as.vector(b1b4.0))
    R = (p00/(1-p00)) * (P2) * ((digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i))^2)
    R2 = t(Ci2) %*% (Ci2 * as.vector(R))
    B22 = b1b4.1 - R2
    #a*b
    b1b5.0 = (trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i)) * (P2)
    b1b5.1 = t(Ci2) %*% (Ci3 * as.vector(b1b5.0))
    O = (p00/(1-p00)) * (P2) * (digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) * (digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i)))
    O2 = t(Ci2) %*% (Ci3 * as.vector(O))
    B23 = b1b5.1 - O2
    #b*b
    b2b6.0 = (-E22 + trigamma(theta1i + theta2i + theta3i) - trigamma(theta2i + theta3i) + trigamma(theta3i)) * (P2)
    b2b6.1 = t(Ci3) %*% (Ci3 * as.vector(b2b6.0))
    Q = (p00/(1-p00)) * (P2) * ((digamma(theta1i+theta3i) + digamma(theta2i+theta3i) - digamma(theta1i+theta2i+theta3i) - digamma(theta3i))^2)
    Q2 = t(Ci3) %*% (Ci3 * as.vector(Q))
    B44 = b2b6.1 - Q2

    F1 = cbind(  A11,    A12,    A13,   A14)
    F2 = cbind(t(A12),   B11,    B12,   B13)
    F3 = cbind(t(A13), t(B12),   B22,   B23)
    F4 = cbind(t(A14), t(B13),  t(B23), B44)

    Fisher=round(rbind(F1, F2, F3, F4), 2)
    Fisher=round(rbind(F1, F2, F3, F4), 2)
    is_positive_definite <- isSymmetric(Fisher) && all(eigen(Fisher)$values > 0)

    if (is_positive_definite) {
      # If Fisher is already positive definite, take its inverse
      Fisher_inv <- solve(Fisher)
    } else {
      # If Fisher is not positive definite, make it positive definite
      Fisher_pos_def <- Matrix::nearPD(Fisher)
      Fisher_pd <- Fisher_pos_def$mat
      Fisher_inv <- solve(Fisher_pd)
    }
    EST<-c(g,b1,b2,b3)
    df = length(EST)
    A1 = EST + 1.96 * sqrt(diag(Fisher_inv))
    A2 = EST - 1.96 * sqrt(diag(Fisher_inv))
    Conf = matrix(c(A2,A1), nrow = (np*2)+2, ncol = 2)
    Dev = A1-A2
    Sig = abs(EST/Dev)
    S_E = sqrt(diag(Fisher_inv)/n)
    z_score <- EST / S_E
    ConfInt = round(cbind(Conf, EST, Sig, S_E, z_score),3)
    colnames(ConfInt) = c("LB", "UB", "Est", "Est / Length", "SE", "z-score")
    return(list(FisherInformation = Fisher, ConfidenceIntervals = ConfInt))
  }
}
