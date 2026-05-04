---
title: "Code to simulate incidence-proxy data"
output:
  rmarkdown::html_vignette:
    keep_md: true
  rmarkdown::github_document:
    toc: true
bibliography: ref.bib
vignette: >
  %\VignetteIndexEntry{simulate-data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




To simulate data, we use a subset of 704 contigous french postal
codes, among which 183 were supposedly covered by a cancer
registry. The shape file is available in the `shap_geo` file. 
Population data by age for these postal codes as well as incidence
rate, se and ppv values by age estimated for Corpus Uteri cancer
in women (CU) with the hospitalisation proxy are given in dataset `dt_cu`. 

The parameters used for the simulation were the same as those used in
the simulation section of the manuscript for the $\rho=0.7$ and $N=20$
scenario. They are given in the following table:


```
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.2.1     ✔ readr     2.2.0
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
#> ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
#> ✔ purrr     1.2.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
#> 
#> Attaching package: 'kableExtra'
#> 
#> 
#> The following object is masked from 'package:dplyr':
#> 
#>     group_rows
```

<table class="table table-striped table-hover" style="color: black; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> param </th>
   <th style="text-align:right;"> value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> $N$ </td>
   <td style="text-align:right;"> 20.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\overline{Se}$ </td>
   <td style="text-align:right;"> 0.80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\overline{PPV}$ </td>
   <td style="text-align:right;"> 0.70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\sigma_c$ </td>
   <td style="text-align:right;"> 0.20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\phi_c$ </td>
   <td style="text-align:right;"> 0.70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\delta_1$ </td>
   <td style="text-align:right;"> 0.57 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\phi_1$ </td>
   <td style="text-align:right;"> 0.70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\delta_2$ </td>
   <td style="text-align:right;"> 0.41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\phi_2$ </td>
   <td style="text-align:right;"> 0.70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\rho$ </td>
   <td style="text-align:right;"> 0.70 </td>
  </tr>
</tbody>
</table>





We first load required R packages, some utility fonctions needed for
simulations and these datasets.


``` r
## 
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(spdep)
#> Loading required package: spData
#> To access larger datasets in this package, install the spDataLarge
#> package with: `install.packages('spDataLarge',
#> repos='https://nowosad.github.io/drat/', type='source')`
library(RSpectra) 
library(Matrix)
#> 
#> Attaching package: 'Matrix'
#> The following objects are masked from 'package:tidyr':
#> 
#>     expand, pack, unpack
library(tidyverse)
library(magrittr)
#> 
#> Attaching package: 'magrittr'
#> The following object is masked from 'package:purrr':
#> 
#>     set_names
#> The following object is masked from 'package:tidyr':
#> 
#>     extract
library(here)
#> here() starts at /home/onyxia/work/bpsp
source(here("R/format_data.R")) ## R utilities to simulate data

## Shapefile
load(file=here("data/shap_geo.rda"))
## Data with pop, muI rate and se and ppv estimated 
load(file=here("data/dt_cu.rda"))
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

The simulation then consists in:

1. Scale the $\mu^I_{ij}$, $Se_{ij}$ and $PPV_{ij}$ levels in `dt_cu` to match the values
for our simulation scenario (i.e., $N$ mean number of total cases in
   geographical areas, $\overline{Se}$ and $\overline{PPV}$);
2. Simulate BYM2 distributed random effects $c_j$,  $u_{1j}$ and
   $u_{2j}$ and update $\mu^I_{ij}$, $Se_{ij}$ and
   $PPV_{ij}$ accordingly;
3. Calculate $\lambda_{ij}^{FN}$, $\lambda_{ij}^{FP}$ and
   $\lambda_{ij}^{TP}$ from $\mu^I_{ij}$, $Se_{ij}$ and $PPV_{ij}$;
4. Generate $FN_{ij}$,  $FP_{ij}$ and $TP_{ij}$ from Poisson
   distribution with mean $\lambda_{ij}^{FN}$, $\lambda_{ij}^{FP}$ and
   $\lambda_{ij}^{TP}$;
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
set.seed(3)

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
#> 1  20.4  23.3   4.33   7.21   16.1 0.789 0.695

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


