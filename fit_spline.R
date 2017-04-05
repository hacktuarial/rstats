fit_spline <- function(data, x, y, df, bs='tp', ...) {
    library(mgcv)
    " Fit a spline like y ~ s(x) using mgcv "
    " Returns: table like |x|x_1|x_2|...|x_k| to merge in " 
    ff <- paste0(y, " ~ s(", x, ", fx=TRUE, bs='", bs, "', k=", df+1, ")")
    m <- gam(as.formula(ff), data=data, family=gaussian, ...)
    # by default, model fits an intercept. drop this
    spline_matrix <- predict(m, type="lpmatrix")[, 2:(df+1)]
    # and add the original x back
    spline_matrix <- cbind(m$model[[x]], spline_matrix)
    colnames(spline_matrix) <- paste0(x, c("", 1:df))
    spline_matrix <- unique(spline_matrix)
    as.data.frame(spline_matrix)
}
