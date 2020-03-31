
# which_constant
.which_constant <- function(x, tol=1e-4){
  which(apply(x, 2, stats::sd) < tol)
}

.which_collinear <- function(x, tol=1e-4){
  cor_m <- stats::cor(x)
  cor_m[upper.tri(cor_m, diag = TRUE)] <- 0
  apply(cor_m, 2, function(.x) any(abs(.x) > (1-tol))) %>% which
}

