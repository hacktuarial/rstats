# calculate bias-corrected and accelerated bootstrap
# confidence intervals
# code written by Prof Efron, Stanford University
# for STATS306A, Winter 2015
# comments by T Sweetser, @hacktuarial
bca <- function (xx, t0, z0, a=0, ASL, 
                 al = c(0.025, 0.05, 0.1, 0.16, 0.84, 0.9, 0.95, 0.975)) {
  ## INPUTS
  # xx = vector of bootstrap estimates
  # t0 = maximum likelihood point estimate
  # z0 = bias correction. calculated automatically if missing
  # a is acceleration. also calculated automatically if value==99
  # not sure what ASL does
  
  ## OUTPUTS
  # 1. Call
  # 2. Confidence intervals
  # 3. theta hat, acceleration, bias correction
  # 4. standard deviation, mean, and # of bootstrap estimates
  # 5. ???
    call <- sys.call()
    B <- length(xx)
    sd.mean.B <- c(sd(xx), mean(xx), B)
    if (missing(z0)) { z0 <- qnorm(mean(xx < t0)) }
    Z <- z0 + qnorm(al)
    if(a==99) {
      xlo <- quantile(xx,.005)
      xhi <- quantile(xx,.995)
      xxx <- pmax(pmin(xx,xhi), xlo)
      xxm <- mean(xxx)  
      a <- mean((xxx-xxm)^3)/(sd(xxx)^3)/6
      # nonparametric skewness. p 5.16 of http://web.stanford.edu/~swager/306a/Winter2015_Stats306A_Notes5.pdf
    }         


    atil <- pnorm(z0 + Z/(1 - a * Z))
    lims <- quantile(xx, atil)
    a0 <- pnorm(z0)
    Zhat <- Z - z0
    Ztil <- qnorm(atil)
    R <- ((Zhat + z0)/(Ztil - z0))^2
    W <- exp(0.5 * (Ztil^2 - Zhat^2)) * R
    if (!missing(z0)) 
        V <- 0
    else V <- -(1 + R) * exp(0.5 * (z0^2 - Zhat^2))
    cv <- a0 * (1 - a0) * V^2 + atil * (1 - atil) * W^2
    cv <- cv + 2 * pmin(a0, atil) * (1 - pmax(a0, atil)) * V * W
    cv <- sqrt(cv/B)/pmin(al, 1 - al)
    lims <- rbind(lims,atil,cv)
    dimnames(lims) <- list(c("bcalims","%iles","CoefVar"),al)
    if (missing(t0)) thetahat <- NA
    else thetahat <- t0
    thet.a.z0 <- c(thetahat, a, z0)
    u1 <- mean(xx)
    u2 <- mean((xx - u1)^2)
    u4 <- mean((xx - u1)^4)
    sdCV <- sqrt((u4/u2^2 - 1)/(4 * length(xx)))
    vl <- list(call, lims, thet.a.z0, sd.mean.B, sdCV)
    if (!missing(ASL)) {
        z1 <- cumsum(.cou(xx, ASL))/length(xx)
        z1 <- qnorm(z1[-(length(z1))]) - z0
        asl <- pnorm(z1/(1 + a * z1) - z0)
        asl <- rbind(ASL, asl)
        vl$asl <- asl
    }
    vl
}
