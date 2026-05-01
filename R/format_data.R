#' inverse logit function
#'
#' @param x 
#'
#' @returns real
#' @export
#'
#' @examples
inv_logit <- function(x) {
  exp(x)/(1+exp(x))
}
#' logit function
#'
#' @param x 
#'
#' @returns real
#' @export
#'
#' @examples
logit <- function(x) {
  log(x/(1-x))
}


##' Scale Q matrix for bym2
##'
##' Take a shape object, and turn it to an adgency matrix,
##' with scaling factor to be used in BYM2 formuation 
##' 
##' @param shape : shapefile
##' @return A named list of class scale_geo with
##' * scales : scaling factors for each connected subgraphs
##' * adj_mat : adjency matrix
##' * dt_comp : tibble of identifiers of geographical area and the subgraph they belong to
##' * Q_gen : matrice de la racine carrée de Q_inv pour générer des EA ICAR
##' 
##' @author Edouard Chatignoux
scale_shap <- function(shape){
  nb<-spdep::poly2nb(shape)
  comps_nb <- spdep::n.comp.nb(nb)
  
  ## Reorder the graph into sub-groups of connected graphs
  dt_comp <- tibble(comp=comps_nb$comp.id,geo=shape$geo)%>%
    group_by(comp)%>%
    mutate(n=n())%>%
    ungroup() %>%
    mutate(id_geo = row_number())
  
  adj_mat <- spdep::nb2mat(nb,style="B", zero.policy = TRUE)%>%
    (Matrix::Matrix) %>%
    as(.,"generalMatrix")
  
  
  sc <- tibble(comp=1:comps_nb$nc)%>%
    mutate(scales=map(comp,function(comp){
      if (sum(dt_comp$comp==comp)>1){
        adj_comp <- adj_mat[dt_comp$comp == comp , dt_comp$comp == comp]
        Q <- Matrix::Diagonal(x=Matrix::rowSums(adj_comp)) - adj_comp %>%
          as(.,"generalMatrix")
        rg<-(ncol(Q)-1)
        if (rg>2)
          QQ<-RSpectra::eigs_sym(Q, k= rg)
        else {
          QQ <- eigen(Q)
          QQ$values <- QQ$values[1:rg]
          QQ$vectors <- QQ$vectors[1:rg,]
        }
        Q_inv <- (QQ$vectors %*% Matrix::Diagonal(x=1/(QQ$values))%*%t(QQ$vectors))
        sc<- Q_inv %>%
          Matrix::diag(.) %>%
          log() %>%
          mean() %>%
          exp()
        Q_gen <- QQ$vectors %*% Matrix::Diagonal(x=1/(sqrt(QQ$values*sc)))
      } else {
        sc = 1
        Q = Matrix::Matrix(0)
        Q_gen = Matrix::Matrix(0)
      }
      tibble(sc=sc, Q_sc=list(Q*sc), Q_gen = list(Q_gen))
    })) %>%
    unnest(scales)
  
  ord_Q<-order(dt_comp %>% arrange(comp,id_geo)%$%id_geo)
  
  ret <- list(scales=sc %>% select(comp,scales=sc),
              dt_comp=dt_comp,
              adj_mat = adj_mat,
              Q_gen = Matrix::bdiag(sc$Q_gen)[ord_Q,]
  )
  class(ret) <- "scale_geo"
  ret
}

##' List of edges
##'
##' List the edges in a scale_geo
##' 
##' @param scales_geo : a list created with [scale_geo]
##' @return A tibble of connected nodes nodes1, nodes2
##' @author Edouard Chatignoux
edges <- function( scales_geo ){
  ## Graph as node pairs
  nod_pairs<-scales_geo$adj_mat%>%
    as.matrix%>%
    as_tibble%>%
    mutate(id=row_number())%>%
    gather(var,val,-id)%>%
    filter(val>0)%>%
    mutate(id2=as.numeric(str_remove(var,"V")))%>%
    arrange(id)%>%
    select(id,id2)%>%
    filter(id<id2)%>%
    rename(nodes1=id,nodes2=id2)
  
  nod_pairs%>%
    left_join(scales_geo$dt_comp%>%rename(nodes1 = id_geo,geo1=geo), by="nodes1")%>%
    left_join(scales_geo$dt_comp%>%rename(nodes2 = id_geo,geo2=geo),by=c("nodes2","comp","n"))%>%
    select(comp,n,nodes1,nodes2,geo1,geo2) 
}


#' Turn geographical neighbourhood structure into stan data for BP-SP usaae
#'
#' @param scales_geo: a scale_geo type object created with [scale_geo]
#' @param sub_ra: character vector containing the names of geographical aeras covered by a registry 
#'
#' @returns
#' @export
#'
#' @examples
prep_geo <- function( scales_geo){
  nodes <- edges(scales_geo)
  n_geo <- scales_geo$dt_comp %>% nrow
  n_comp <- scales_geo$scales$comp %>% max
  
  comp_idx <-
    scales_geo$dt_comp %>% 
    arrange(comp)%$%id_geo
  comp_size <- scales_geo$dt_comp %>% 
    group_by(comp) %>%
    summarise(n=n())%$%n
  
  ret<-list(
    n_geo = n_geo,
    n_comp = n_comp,
    comp_size = array(comp_size, dim = n_comp), 
    n_edges = nrow(nodes), 
    node1 = nodes$nodes1, 
    node2 = nodes$nodes2, 
    comp_idx = array(comp_idx, dim = n_geo),
    scales = scales_geo$scales$scales %>% as.array)
  class(ret) <- "geo_stan"
  ret

}

#' Turn smooth specification to equivalent random spe, for stan use
#'
#' @param smooth : an mgcv like spline specification
#' @param data : the data to evaluate the spline
#' @param newdata : if spline has to be predicted on new data
#'
#' @returns The matrix of spline in re format
#' @export
#'
#' @examples
#' smooth2stan(s(x,bs="tp"), data=tibble(x=1:100))
smooth2stan <- function(smooth,data,newdata = NULL){
  sm <- mgcv:::smoothCon(smooth,
                         absorb.cons=T,
                         modCon = 3,
                         data=data)[[1]]
  re <- mgcv:::smooth2random(sm,"",2)
  if (!is.null(newdata)) {
    X <- PredictMat(sm,newdata)   ## get prediction matrix for new data
    ## transform to r.e. parameterization
    if (!is.null(re$trans.U)) X <- X%*%re$trans.U
    X <- t(t(X)*re$trans.D)
    ## re-order columns according to random effect re-ordering...
    X[,re$rind] <- X[,re$pen.ind!=0] 
    ## re-order penalization index in same way  
    pen.ind <- re$pen.ind; pen.ind[re$rind] <- pen.ind[pen.ind>0]
    ## start return object...
    r <- list(rand=list(),Xf=X[,which(re$pen.ind==0),drop=FALSE])
    for (i in 1:length(re$rand)) { ## loop over random effect matrices
      r$rand[[i]] <- X[,which(pen.ind==i),drop=FALSE]
      attr(r$rand[[i]],"s.label") <- attr(re$rand[[i]],"s.label")
    }
    names(r$rand) <- names(re$rand)
    re<-r
  }
  cbind(re$Xf,re$rand[[1]]) 
}

##' Extract components from the formula 
##'
##' @title  get_comp
##' @param form : a formulae
##' @return A tibble that contains the compoennet of the formulae.
##' @details This is an internal function intended to be used in prep_data
##' @author Edouard Chatignoux
get_comp <- function(form){
  
  rhs<-labels(terms(form))%>%as_tibble() %>%
    rename(call=value)
  rhs %<>% 
    rowwise() %>%
    mutate(prox = str_extract(call,"[a-z]+"),
           for_all = !(prox %in% c("i","se","ppv")),
           effects = ifelse(for_all, 
                            call,
                            str_split_1(call,paste0(prox,"\\("))[2] %>%
                              str_remove("\\)")),
           effects=str_remove_all(effects," ") %>%
             str_split("\\+")
           ) %>%
    unnest(effects) %>%
    rowwise() %>%
    mutate(type =
             case_when(
               str_detect(effects,"^s\\(") & !(str_detect(effects,"bym|icar|iid"))~"spline",
               str_detect(effects,"bym|icar|iid")~"re",
               TRUE ~ "fixe"
             )) %>%
    mutate(value =
             case_when(
               str_detect(effects,"bym")~3,
               str_detect(effects,"iid")~1,
               str_detect(effects,"icar")~2,
               type=="spline" & !str_detect(effects,"fx=TRUE")~ 1,
               TRUE ~ 0
             )) 
  rhs %<>% 
    mutate(prox=ifelse(for_all,list(c("i","se","ppv")),list(prox))) %>%
    unnest(prox) %>%
    rowwise() %>%
    mutate(cov=ifelse(str_detect(effects,"\\("),
             str_extract(str_split_1(effects,"\\("),"[a-z]+")[2],
             effects
             )) %>%
    mutate(prox=factor(prox,levels=c("i","se","ppv"))) %>%
    arrange(prox) 
  
  rhs %>%
    arrange(prox,cov,for_all) %>%
    group_by(prox,cov) %>%
    slice(1) %>%
    mutate(effects = str_remove(effects,",fx=TRUE") %>%
             str_replace("\\)\\)","\\)"))%>%
    select(-call,-for_all)

}

#' Covert formula for covariates effects to stan compatible data
#'
#' @param form: formula for effects shared by I, Se and PPV
#' @param I: formula for effects specific to I
#' @param se: formula for effects specific to Se 
#' @param ppv: formula for effects specific to PPV 
#' @param dt_cov: tibble containing covariates values for each geographical areas
#' @param newdata: if prediction on new data
#'
#' @returns list
#' @export
#'
#' @examples
prep_cov <- function(
    form = ~s(chom,k=5)+s(apl,k=5),
    dt_cov=dt_cov,
    newdata = NULL
){

  if (is.null(newdata))
    newdata <- dt_cov
  cov <- get_comp(form) %>% filter(type!="re")
  ## Covars
  W <- cov %>%
    ungroup() %>%
    select(effects) %>%
    unique() %>%
    mutate( W = map(effects,~{
      if (str_detect(.x,"s\\(|t2\\(")) {
        W <- smooth2stan(eval(parse(text=.x)),
                          data=dt_cov,
                         newdata=newdata
                         )
      } else {
        W<-model.matrix(as.formula(paste0("~-1+",.x)),data=newdata)
      }
    })) 
  W  %<>% arrange(effects)
  
  has_w <- cov %>% ungroup() %>% select(prox,effects,value) %>% 
    pivot_wider(names_from = effects, values_from = value) %>%
    mutate(across(-prox,~as.numeric(!is.na(.x))))

  is_pen_w <- cov %>%
    ungroup() %>%
    select(effects,value) %>%
    unique() %$% value %>% as.array()

  n_w <- length(unique(cov$effects))
  has_w <- has_w %>% select(-prox) %>% as.matrix %>% t()
  d_w <- W  %>% mutate(ncol=map_dbl(W,ncol))%$%ncol %>% as.array()
  W <- W  %$% bind_cols(W)
  if ( ncol(W)==0 ) {
    n_w=0
    W<-matrix(ncol=0,nrow=nrow(cov)) 
  }
  list(dt_cov=dt_cov,
       d_W=ncol(W),
       W=W,
       n_w=n_w,
       has_w=has_w,
       d_w=d_w,
       is_pen_w =is_pen_w
  )
}


##' Prepare data to be used in stan fit
##'
##' @param geo : geographical info (sf, scale_geo or geo_stan object)
##' @param dt : data frame with geo 
##' @param form : formula
##' @return A list suited for stan sampling
##' @author Edouard Chatignoux
prep_data <- function(dt,
                      form_id = c(k,ald) ~ offset(log(pa)) + s(age,bs="cr",k=5,fx=TRUE),
                      form_geo= ~ bym(geo),
                      geo = shape_geo
                      ){
  
  ## Neigbour structure 
  if (!is.null(geo)){
    if (class(geo)[1]!="geo_stan") {
      ## Geography for stan use
      geo_stan <- prep_geo(scale_shap(geo))
    } else {
      geo_stan <- geo
    }
  } else if (class(dt)[[1]]=="sf"){
    geo_stan <- dt %>% 
      select(geo,geometry) %>% 
      unique() %>%
      arrange(geo)
    geo_stan <- prep_geo(scale_shap(geo_stan))
  } else {
    warning("No geo supplied")
    geo_stan<-list()
  }
    
  ## Names of counting vars and incidence and proxy labels
  cts <- rlang::f_lhs(form_id) %>% all.vars()
  inc <- cts[1]
  prox <- cts[2]

  ## Decode formulaes
  cov_id <- get_comp(form_id) %>%
    filter(prox=="i") %>%
    ungroup() %>%
    select(-prox)
  cov_geo <- get_comp(form_geo) 
  

  ## Select data
  stan_dt <- dt %>% ungroup()
  
  ## Assign offsets
  off<-terms(form_id) %>% attr(.,"offset")
  if (is.null(off)){
    stan_dt[,"off"]<-1
  }  else{
    name_off <- all.vars(form_id[-2])[off-1]
    off<-attributes(terms(form_id))$variable %>% as.character %>% .[off+1]
    off <- model.frame( as.formula(str_glue("~{f}",f=off)), data=stan_dt)
    stan_dt[,"off"]<-off
    stan_dt %<>% select(-any_of(name_off))
  }

  
  ## Select working variables
  stan_dt %<>%
    select(geo,RA,age=!!cov_id$cov,!!cov_geo$cov,inc=!!inc,prox=!!prox,off)

  ## Geos in RA
  geo_ra <- stan_dt %>% filter(!is.na(inc)) %$% geo %>% unique()
  stan_dt %<>%
    mutate( RA = geo %in% geo_ra) %>%
    mutate(inc=replace_na(inc,999999)) # stan doesn't like NA
    
  ## Ids for geo and age groups
  stan_dt %<>% 
    arrange(age,geo) %>%
    mutate(id_geo = as.numeric(as.factor(geo)),
           id_age = as.numeric(as.factor(age)),
           id_data = row_number()
    ) %>%
    ungroup() 
  
  ## Age effect
  Z <- smooth2stan(eval(parse(text=cov_id$effects)),
                   data=stan_dt%>% select(age) %>% unique)

  ## Covariates
  cov_stan <- prep_cov(form = form_geo,
                       dt_cov=stan_dt %>% ungroup() %>% select(geo,!!cov_geo$cov) %>% unique())
    
  
  decod <- list(age=stan_dt %>% select(id_age,age) %>% unique,
                geo=stan_dt %>% select(id_geo,geo,RA) %>% unique ,
                geo_ra=stan_dt %>% filter(RA) %>%
                  select(geo,RA) %>% unique %>%
                  mutate(id_geo=as.numeric(as.factor(geo))))
                
  stan_dt <- 
    list(prior_only = 0,
         n_obs = nrow(stan_dt),
         inc = stan_dt %>% select(inc,prox) %>% as.matrix,
         off = stan_dt$off,
         n_age = max(stan_dt$id_age),
         id_age=stan_dt$id_age,
         d_z=ncol(Z),
         Z=Z,
         id_geo=stan_dt$id_geo,
         n_geo_ra =  length(geo_ra),
         n_geo_ora = max(stan_dt$id_geo) -length(geo_ra),
         geo_ra = stan_dt %>% filter(RA) %$% id_geo %>% unique,
         geo_ora = stan_dt %>% filter(!RA) %$% id_geo %>% unique,
         id_geo_ra=stan_dt%>% filter(RA) %>% mutate(id=as.numeric(as.factor(id_geo))) %$% id,
         id_geo_ora=stan_dt%>% filter(!RA) %>% mutate(id=as.numeric(as.factor(id_geo))) %$% id,
         n_obs_ra = stan_dt %>% filter(RA) %>% nrow,
         n_obs_ora = stan_dt %>% filter(!RA) %>% nrow,
         id_obs_ra= stan_dt %>% filter(RA)%$%id_data, 
         id_obs_ora= stan_dt %>% filter(!RA)%$%id_data,
         re_typ=cov_geo %>% filter(cov=="geo")%$%value %>% as.array(),
         ra_only = 0,
         is_pen = cov_id$value,
         is_cor = 1,
         p_ru = c(1.5,1.5) %>% as.array()
    )
  
  tibble(
    dt_fit  = list(c(stan_dt,
                  geo_stan,
                  cov_stan[-1])),
    dt_post = list(tibble(decod=list(decod),
               dt_cov=list(cov_stan$dt_cov))),
    form_id=list(form_id),
    form_geo=list(form_geo))
         
}


