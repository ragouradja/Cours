```{r}
# 1 --> Normalisé par D
estim_pi <- function(n, r) {
  theta = runif(n,0,pi)
  x = runif(n,0,1)
  a = x + cos(theta)* r
  ns = sum(a > 1 | a < 0)
  ps = ns/n
  estimateur_pi = 2*r/ps
  return(estimateur_pi)
}




```

```{r}
# L / D = 2 / 5 = 0.4
estim_pi(10000,0.3)
```


```{r}
require(MASS)
```

```{r}
mse_pi = function(n,r,nrep){
  valeur_pi = c()
  for (run in 1:nrep) {
  valeur_pi = c(valeur_pi, estim_pi(n,r))
  }
  MSE = mean((valeur_pi - pi)**2)
  truehist(valeur_pi)
  abline(v = pi, col = "red")
  print(MSE)
}
```

```{r}
mse_pi(1000,0.3,100000)
```
```{r}
cercle = function(n){
  r = 1
  x = runif(n,-r,r)
  y = runif(n,-r,r)
  points = 0
  d = sqrt((x**2 + y**2))
  list_points = ifelse(d <= r, 1, 0)
  points = sum(list_points)
  plot(x,y)
  print(points/n*4)
}
```

```{r}
cercle(1000)
pi
```

```{r}
pn = rnorm(10000,0,1)
truehist(pn)
mean(pn)
var(pn)

p = rnorm(10000,0,1)**2
truehist(p)
val_x = seq(from = 0, to = 7, len = 10000)
val_y = dchisq(val_x, df= 2)
lines(val_x,val_y, col = "red", lwd = 2)
```
```{r}
quantile(pn,probs = c(0.500112))
```

