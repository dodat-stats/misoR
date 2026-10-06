# Joint population posterior and numerical diagnostics. Gamma uses shape/rate.
.miso_quadrature_control <- function(control=list()) {
  defaults <- list(tol=1e-8,max_order=512L,tail_drop=40,tail_expand=10)
  .miso_stop(is.list(control) && (length(control)==0L ||
    (!is.null(names(control)) && all(nzchar(names(control))) &&
     !anyDuplicated(names(control)) && all(names(control) %in% names(defaults)))),
    "Unknown or unnamed quadrature control.")
  defaults[names(control)] <- control
  .miso_stop(all(vapply(defaults,function(x) is.numeric(x) && length(x)==1L &&
    is.finite(x) && x>0,TRUE)),"Quadrature controls must be positive finite scalars.")
  .miso_stop(defaults$max_order %in% 2^(6:11),"max_order must be a power of two from 64 to 2048.")
  defaults
}

#' Calibrate the population prior from observed row totals
#' @param Y Dense or sparse nonnegative observations by features.
#' @param D Number of slots per motif.
#' @param n Optional positive exposures; use NULL for ordinary MiSO.
#' @param alpha_mean Prior mean of the population shape, default 1.
#' @return A prior specification with named shape/rate pairs `alpha` and `beta`,
#'   and calibration diagnostics. Both shapes are 2. Its expected slot loading
#'   is mean(rowSums(Y)/n)/D (n=1 for ordinary MiSO). The calibration rule stays
#'   fixed across D; numerical rates stay fixed across S and starts.
#' @export
miso_population_prior <- function(Y,D,n=NULL,alpha_mean=1) {
  Y <- .miso_prepare_counts(Y); .miso_positive_integer(D,"D")
  .miso_stop(is.numeric(alpha_mean) && length(alpha_mean)==1L &&
    is.finite(alpha_mean) && alpha_mean>0,"alpha_mean must be positive.")
  .miso_check_exposure(n,nrow(Y))
  total <- mean(.miso_row_sums(Y)/(if(is.null(n)) 1 else n))
  .miso_stop(total>0,"Row-total calibration requires positive total data mass; supply an explicit population_prior.")
  list(alpha=c(shape=2,rate=2/alpha_mean),beta=c(shape=2,rate=total/(D*alpha_mean)),
       calibration=list(method="row_totals",mean_total=total,mean_slot=total/D,
                        alpha_mean=alpha_mean,D=D,units=if(is.null(n)) "counts" else "per_unit_n"))
}

.miso_check_exposure <- function(n,N) {
  if(!is.null(n)) .miso_stop(is.numeric(n) && is.null(dim(n)) && length(n)==N &&
    all(is.finite(n) & n>0),"n must contain one positive finite exposure per observation.")
}

.miso_resolve_prior <- function(prior,Y,D,n) {
  if(is.null(prior) || identical(prior,"row_totals")) prior <- miso_population_prior(Y,D,n)
  .miso_stop(is.list(prior) && all(c("alpha","beta") %in% names(prior)),
    "population_prior must contain alpha and beta shape/rate pairs, or be 'row_totals'.")
  for (name in c("alpha","beta")) {
    x <- prior[[name]]
    .miso_stop(is.numeric(x) && length(x)==2L && setequal(names(x),c("shape","rate")) &&
      all(is.finite(x) & x>0),"Prior pairs need positive named shape and rate.")
    prior[[name]] <- x[c("shape","rate")]
  }
  .miso_stop(prior$beta["shape"]>1,"The beta prior shape must exceed 1 for a finite population loading mean.")
  prior
}
.miso_prior_vector <- function(prior) unname(c(prior$alpha,prior$beta))

.miso_population_update <- function(a,b,omega,prior,control) {
  N <- dim(a)[1]; S <- dim(a)[2]; D <- dim(a)[3]
  blocks <- vector("list",S*D)
  for(s in seq_len(S)) {
    rate <- .miso_rate_matrix(b,s,N,D)
    for(d in seq_len(D)) {
      blocks[[s+(d-1L)*S]] <- .miso_gamma_pair(sum(omega[,s]),
        sum(omega[,s]*a[,s,d]/rate[,d]),
        sum(omega[,s]*(digamma(a[,s,d])-log(rate[,d]))),.miso_prior_vector(prior),control)
    }
  }
  blocks
}
.miso_population_matrix <- function(blocks,field,S,D)
  matrix(vapply(blocks,function(h) unname(h[[field]]),0.0),S,D)
.miso_quadrature_diagnostics <- function(blocks,S,D) {
  data.frame(motif=rep(seq_len(S),D),slot=rep(seq_len(D),each=S),
    converged=vapply(blocks,function(h) h$converged,TRUE),
    order=vapply(blocks,function(h) as.integer(h$quadrature_order),0L),
    refinement_error=vapply(blocks,function(h) h$quadrature_error,0.0),
    tail_check_error=vapply(blocks,function(h) h$tail_check_error,0.0),
    log_alpha_lower=vapply(blocks,function(h) h$bounds[1],0.0),
    log_alpha_upper=vapply(blocks,function(h) h$bounds[2],0.0),
    mode_log_alpha=vapply(blocks,function(h) h$mode_log_alpha,0.0))
}

.miso_hyper_rule <- local({
  cache <- new.env(parent=emptyenv())
  function(n) {
    key <- as.character(n)
    if(!exists(key,cache,inherits=FALSE)) {
      j <- seq_len(n-1L); b <- j/sqrt(4*j*j-1)
      J <- matrix(0,n,n); J[cbind(j,j+1L)] <- b; J[cbind(j+1L,j)] <- b
      ev <- eigen(J,symmetric=TRUE); o <- order(ev$values)
      assign(key,list(x=ev$values[o],w=2*ev$vectors[1,o]^2),cache)
    }
    get(key,cache,inherits=FALSE)
  }
})

# Optimal joint block under independent normalized Gamma hyperpriors.
# n=sum omega, T=sum omega E(lambda), R=sum omega E(log lambda).
# beta | alpha ~ Gamma(C+n*alpha, Bbeta+T); integrate beta analytically.
.miso_gamma_pair <- function(n,T,R,prior,control=list()) {
  control <- .miso_quadrature_control(control)
  tol <- control$tol
  if(length(n)!=1L || length(T)!=1L || length(R)!=1L ||
     !all(is.finite(c(n,T,R))) || n<0 || T<0 ||
     length(prior)!=4L || any(!is.finite(prior) | prior<=0) || prior[3]<=1)
    stop("Invalid joint Gamma block arguments.",call.=FALSE)
  A <- unname(prior[1]); B <- unname(prior[2])
  C <- unname(prior[3]); Db <- unname(prior[4]); v <- Db+T
  # The omitted -R is independent of alpha; restore it in logZ below.
  logkernel <- function(t) {
    a <- exp(t)
    A*t+(R-B)*a-n*lgamma(a)+lgamma(C+n*a)-(C+n*a)*log(v)
  }
  lo <- -25; hi <- 20
  mode <- stats::optimize(logkernel,c(lo,hi),maximum=TRUE,tol=1e-9)$maximum
  if(mode-lo < .1 || hi-mode < .1) stop("Gamma block mode reached search boundary.",call.=FALSE)
  peak <- logkernel(mode)
  interval <- function(drop) {
    tailfn <- function(t) logkernel(t)-peak+drop
    left <- mode-1; right <- mode+1
    for (j in seq_len(100L)) { if(tailfn(left)<=0) break; left <- left-2 }
    for (j in seq_len(100L)) { if(tailfn(right)<=0) break; right <- right+2 }
    c(stats::uniroot(tailfn,c(left,mode),tol=1e-10)$root,
      stats::uniroot(tailfn,c(mode,right),tol=1e-10)$root)
  }
  base_bounds <- interval(control$tail_drop)
  bounds <- interval(control$tail_drop+control$tail_expand)
  quadrature <- function(order, bounds) {
    rule <- .miso_hyper_rule(order)
    left <- bounds[1]; right <- bounds[2]
    t <- (left+right)/2+(right-left)*rule$x/2
    a <- exp(t); u <- C+n*a
    raw <- rule$w*(right-left)/2*exp(logkernel(t)-peak)
    z <- sum(raw); w <- raw/z
    elogb <- digamma(u)-log(v)
    fun <- cbind(mean_alpha=a,mean_beta=u/v,mean_log_alpha=t,
      mean_log_beta=elogb,mean_lgamma_alpha=lgamma(a),
      mean_alpha_log_beta=a*elogb,mean_alpha_over_beta=a*v/(u-1),
      second_alpha=a*a,second_beta=u*(u+1)/v^2,mean_alpha_beta=a*u/v)
    moments <- colSums(fun*w)
    logZ <- peak+log(z)+A*log(B)-lgamma(A)+C*log(Db)-lgamma(C)-R
    list(moments=moments,logZ=logZ)
  }
  discrepancy <- function(x,y) max(abs(x$moments-y$moments)/(1+abs(x$moments)),
                                    abs(x$logZ-y$logZ))
  prev <- quadrature(32L,bounds); error <- tail_error <- Inf
  for(order in 2^(6:log2(control$max_order))) {
    cur <- quadrature(order,bounds)
    base <- quadrature(order,base_bounds)
    error <- discrepancy(cur,prev)
    tail_error <- discrepancy(cur,base)
    if(error<tol && tail_error<tol) break
    prev <- cur
  }
  if(!is.finite(error) || !is.finite(tail_error) || max(error,tail_error)>=tol)
    stop("Population quadrature did not converge; increase max_order or tail_drop.",call.=FALSE)
  m <- cur$moments
  KL <- n*(m["mean_alpha_log_beta"]-m["mean_lgamma_alpha"]) +
    (m["mean_alpha"]-1)*R-m["mean_beta"]*T-cur$logZ
  if(KL < -1e-6) stop("Negative joint Gamma KL: numerical integration failed.",call.=FALSE)
  result <- as.list(m)
  result$KL <- max(0,unname(KL));result$logZ <- cur$logZ
  result$var_alpha <- max(0,unname(m["second_alpha"]-m["mean_alpha"]^2))
  result$var_beta <- max(0,unname(m["second_beta"]-m["mean_beta"]^2))
  result$cov_alpha_beta <- unname(m["mean_alpha_beta"]-m["mean_alpha"]*m["mean_beta"])
  result$n <- n;result$T <- T;result$R <- R;result$prior <- prior
  result$mode_log_alpha <- mode;result$bounds <- bounds
  result$quadrature_order <- order;result$quadrature_error <- error
  result$tail_check_error <- tail_error; result$converged <- TRUE
  result
}

