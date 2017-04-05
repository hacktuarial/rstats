c_test <- function(y, p, g=10) {
    "Hosmer-Lemeshow test for logistic regression goodness-of-fit"
    " y: vector of binary responses
      p: fitted probabilities, e.g. from glm
      g: number of groups. Default: 10 "
    # check arguments
    stopifnot(length(unique(y)) == 2)
    stopifnot(min(p) >= 0)
    stopifnot(max(p) <= 1)
    stopifnot(g > 2)
    stopifnot(is.integer(g))
    G <- gtools::quantcut(p, q=g, labels=1:g)
    p_bar <- purrr::map_dbl(split(p, G), mean)
    obs <- purrr::map_int(split(y, G), sum)
    n <- table(G)
    expected <- p_bar * n
    C <- sum((obs - expected)**2 / (n * p_bar * (1 - p_bar))) 
    p_value <- pchisq(C, df=g-2, lower.tail=FALSE)
    list(C=C, p=p_value)
}
