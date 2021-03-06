<!-- README.md is generated from README.Rmd. Please edit that file -->
``` r
joint_drift_model <- function(mu, Sigma) {
  structure(list(mu = mu, Sigma = Sigma),
            class = "drift_model")
}

condition <- function(model, index, freq) {
  mu_1 <- model$mu[-index]
  mu_2 <- model$mu[index]
  Sigma_11 <- model$Sigma[-index, -index]
  Sigma_12 <- model$Sigma[-index, index]
  Sigma_21 <- model$Sigma[index, -index]
  Sigma_22 <- model$Sigma[index, index]
  
  X <- Sigma_12 %*% solve(Sigma_22)
  new_mu <- mu_1 + X %*% (freq - mu_2)
  new_Sigma <- Sigma_11 - X %*% Sigma_21
  
  joint_drift_model(new_mu, new_Sigma)
}

plot.drift_model <- function(x, ..., indices = c(1, 2)) {
  i <- indices[1]
  j <- indices[2]
  mu <- x$mu[c(i,j)]
  Sigma <- x$Sigma[c(i,j), c(i,j)]
  
  range1 <- range2 <- seq(0, 1, length.out = 200)
  data.grid <- expand.grid(f1 = range1, f2 = range2)
  q.samp <- cbind(data.grid, 
                  prob = mvtnorm::dmvnorm(data.grid, mean = mu, sigma = Sigma))

  ggplot(q.samp, aes(x = f1, y = f2, z = prob)) + 
    geom_contour() +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), ratio = 1)
}
```

``` r
test_mu <- c(0.5, 0.2, 0.6)
test_Sigma <- matrix(c(
  0.1,  0.06, 0.1,
  0.06, 0.1,  0.0,
  0.1,  0.0,  0.3
  ), byrow = TRUE, nrow = 3)
model <- joint_drift_model(test_mu, test_Sigma)
plot(model, xlab = "Pop 1", ylab = "Pop 2")
```

![](README-unnamed-chunk-2-1.png)

``` r
plot(model, indices = c(1,3), xlab = "Pop 1", ylab = "Pop 3")
```

![](README-unnamed-chunk-2-2.png)

``` r
plot(model, indices = c(2,3), xlab = "Pop 2", ylab = "Pop 3")
```

![](README-unnamed-chunk-2-3.png)

``` r
plot(condition(model, index = 1, freq = 0.1))
```

![](README-unnamed-chunk-3-1.png)

``` r
plot(condition(model, index = 1, freq = 0.4))
```

![](README-unnamed-chunk-3-2.png)

``` r
plot(condition(model, index = 1, freq = 0.6))
```

![](README-unnamed-chunk-3-3.png)

``` r
plot(condition(model, index = 1, freq = 0.8))
```

![](README-unnamed-chunk-3-4.png)
