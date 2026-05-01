---
title: "Code to simulate incidence-proxy data"
output: rmarkdown::html_vignette
bibliography: ref.bib
vignette: >
  %\VignetteIndexEntry{simulate-data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




``` r
library(sf)
library(tidyverse)
library(magrittr)
source("./R/format_data.R") ## R utilities to format data for BP-SP evaluation
#> Warning in file(filename, "r", encoding = encoding): impossible d'ouvrir le fichier './R/format_data.R' : No
#> such file or directory
#> Error in file(filename, "r", encoding = encoding): impossible d'ouvrir la connexion
```

To simulate data, we use a subset of 704 contigous french postal
codes, among which 183 were supposedly covered by a cancer
registry. The shape file is available in the `shap_geo` file. 
Population data by age for these postal codes as well as incidence
rate, se and ppv values by age estimated for Corpus Uteri cancer
in women (CU) with the hospitalisation proxy are given in dataset `dt_cu`. 

We first load these datasets.


``` r
## Shapefile
load(file="./data/shap_geo.rda")
#> Warning in readChar(con, 5L, useBytes = TRUE): impossible d'ouvrir le fichier compressé
#> './data/shap_geo.rda', cause probable : 'No such file or directory'
#> Error in readChar(con, 5L, useBytes = TRUE): impossible d'ouvrir la connexion
## Data with pop, muI rate and se and ppv estimated 
load(file="./data/dt_cu.rda")
#> Warning in readChar(con, 5L, useBytes = TRUE): impossible d'ouvrir le fichier compressé './data/dt_cu.rda',
#> cause probable : 'No such file or directory'
#> Error in readChar(con, 5L, useBytes = TRUE): impossible d'ouvrir la connexion
```

To generate `ICAR` random effects, we use the relation found for example in @Riebler2016:

$$\boldsymbol{u} \sim N(0,\boldsymbol{Q^{-}}))$$

where $\boldsymbol{Q^{-}}$ is a generalized inverse  of the
neighbouring matrix $\boldsymbol{Q}$, defined as:

$$Q_{ij} =
 \begin{cases}
   \text{Number of neigbours for $i$ if } i=j\\
   -1 \text{ if $i$ and $j$ share a common border}\\
   0 \text{ otherwise}
 \end{cases}.
 $$.
 
Using a SVD decomposition, $\boldsymbol{Q^{-}}$ is decomposed as:

$$\boldsymbol{Q^{-}}=\boldsymbol{Q_g}\boldsymbol{Q_g^\top}$$

with $\boldsymbol{Q_g}$ of dimension $(J,J-1)$, where $J$ is the total
number of geographical areas.

$\boldsymbol{u}$ is then generated as:

$$\boldsymbol{u} = \boldsymbol{Q_g} \boldsymbol{b}$$

where $\boldsymbol{b}$ is as standardized gaussian random variable of
size $J-1$. 

The $\boldsymbol{Q_g}$ matrix is made available trough the function `scale_shap`.


``` r

cp <- scale_shap(shap_geo)

```

The rest of the simulation consists in:

1. Scaling the LM $\mu^I_{ij}$, $Se_{ij}$ and $PPV_{ij}$ levels to match the values
for our simulation scenario (N=20 mean number of total cases in
   geographical areas, $\overline{Se}=0.8$ and $\overline{PPV}=0.7$;
2. Simulate BYM2 distributed random effects $c_j$,  $u_{1j}$ and
   $u_{2j}$ and integrate them to $mu^I_{ij}$, $Se_{ij}$ and
   $PPV_{ij}$;
3. Calculate $\lambda_{ij}^FN$, $\lambda_{ij}^FP$ and
   $\lambda_{ij}^TP$ from $c_j$,  $u_{1j}$ and
   $u_{2j}$ and integrate them to $mu^I_{ij}$, $Se_{ij}$;
4. Generate $FN_{ij}$,  $FP_{ij}$ and $TP_{ij}$ from Poisson
   distribution.
5. Calculate $I_{ij}=FN_{ij}+TP_{ij}$ and $P_{ij}=FP_{ij}+TP_{ij}$.


Details of the code is given bellow:


``` r

###' params ===============
## Simulation parameters
N<-20
se=.8
ppv=.7
s_c =.2
##' Set s_u1 and s_u2 so that overall error in proxy s_e
##' is .1
r_u = 0.7
s_u1 = .57
s_u2 = .41

## Check : 1st order dl of s_e
sqrt((1-se)^2*s_u1^2+(1-ppv)^2*s_u2^2-2*r_u*((1-se)*s_u1*(1-ppv)*s_u2))
#> [1] 0.09216398
## Prop of spatial variance
phi_c=phi_1=phi_2<-.7

n <- nrow(shap_geo)
n_b <- n-1

RNGkind("Wich")
set.seed(2)

## Generate icar components
c_icar = cp$Q_gen %*% rnorm(n_b) %>% as.numeric
u1_icar = cp$Q_gen %*% rnorm(n_b) %>% as.numeric
u2_icar = cp$Q_gen %*% rnorm(n_b) %>% as.numeric
## Generate iid components
c_iid = rnorm(n)
u1_iid = rnorm(n)
u2_iid = rnorm(n)
## Add correlation between u_i's
u2_icar <- r_u * u1_icar + sqrt(1-r_u^2) * u2_icar
u2_iid = r_u * u1_iid + sqrt(1-r_u^2) * u2_iid

## Store bym2 effects in a tibble
dt_re <- shap_geo %>%
  as_tibble() %>%
  select(geo) %>% 
  mutate( c = sqrt(phi_c) * c_icar + sqrt(1-phi_c) * c_iid,
          u1 = sqrt(phi_1) * u1_icar + sqrt(1-phi_1) * u1_iid,
          u2 = sqrt(phi_2) * u2_icar + sqrt(1-phi_2) * u2_iid
  )

dt_sim <- dt_cu %>%
  mutate(muP=se/ppv*muI)
## Scale muI to have a mean N cases per geo
dt_sim %<>%
  group_by(geo) %>%
  mutate(w=sum(muI*py)) %>%
  ungroup() %>%
  mutate(w=N/mean(w)) %>%
  mutate(muI=muI*w) %>%
  select(-w)

dt_sim %>% 
  group_by(geo) %>%
  mutate(w=sum(muI*py)) %$%mean(w)
#> [1] 20

## Scale Se to have mean Se per geo
dt_sim %<>%
  group_by(age) %>%
  mutate(w=sum(muI*py)) %>%
  ungroup() %>%
  mutate(w=w/sum(muI*py)) %>%
  group_by(geo) %>%
  mutate(se0=sum(se*w)) %>%
  mutate(se=se/se0*!!se) %>%
  select(-w,-se0)

## Same for ppv
dt_sim %<>%
  group_by(age) %>%
  mutate(w=sum(muP*py)) %>%
  ungroup() %>%
  mutate(w=w/sum(muP*py)) %>%
  group_by(geo) %>%
  mutate(ppv0=sum(ppv*w)) %>%
  mutate(ppv=ppv/ppv0*!!ppv) %>%
  select(-w,-ppv0)

## i, se, vpp, lbd_fn, lbd_fp and lbd_tp 
dt_sim %<>%
  left_join(dt_re, by = "geo") %>%
  group_by(RA) %>%
  mutate(across(c(u1,u2),~.x-mean(.x))) %>%
  ungroup() %>%
  mutate(
    muI = muI*exp(s_c*c),
    se = inv_logit( logit(se) + s_u1*u1),
    ppv = inv_logit( logit(ppv) + s_u2*u2),
    lbd_fn = muI * (1-se),
    lbd_fp = muI * (1-ppv) * se/ppv,
    lbd_tp = muI * se,
    muP=lbd_fp+lbd_tp
  )

dt_sim %>% 
  group_by(geo) %>%
  mutate(across(c(muI,muP,lbd_fn,lbd_fp,lbd_tp),~sum(.x*py)),
         se = lbd_tp/muI,ppv=lbd_tp/muP) %>%
  ungroup() %>%
  summarise(across(c(muI,muP,lbd_fn,lbd_fp,lbd_tp,se,ppv),mean))
#> # A tibble: 1 × 7
#>     muI   muP lbd_fn lbd_fp lbd_tp    se   ppv
#>   <dbl> <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1  20.8  23.9   4.21   7.31   16.6 0.782 0.695

# Counts TP, FN, FP and I and P
dt_sim %<>%
  rowwise() %>%
  mutate(
    fn=rpois(n=1, lambda = py*lbd_fn),
    fp=rpois(n=1, lambda = py*lbd_fp),
    tp=rpois(n=1, lambda = py*lbd_tp),
    I=fn+tp,
    P=fp+tp
  ) %>% 
  ungroup()

# Set I to NA outside registry area
dt_sim %<>%
  mutate(I=ifelse(RA,I,NA))

```


```
#> Warning in gzfile(file, "wb"): impossible d'ouvrir le fichier compressé './data/dt_sim.rda', cause probable
#> : 'No such file or directory'
#> Error in gzfile(file, "wb"): impossible d'ouvrir la connexion
```
