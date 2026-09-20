## Bound identities and streamed sufficient statistics against independent sums.
local({
  close = function(x, y, tolerance = 1e-10) {
    result = all.equal(x, y, tolerance = tolerance)
    if (!isTRUE(result)) stop(paste(result, collapse = "\n"))
  }
  scores = rbind(c(1000, 999, -1000), c(-1000,-1001,-2000), c(0,-2,-1000))
  reference = apply(scores, 1, function(x) max(x) + log(sum(exp(x-max(x)))))
  close(.miso_log_sum_exp_rows(scores), reference)
  close(.miso_log_sum_exp_rows(matrix(c(-Inf,Inf), 2, 1)), c(-Inf,Inf))
  stopifnot(identical(.miso_log_sum_exp_rows(matrix(numeric(), 0, 3)), numeric()))

  set.seed(1809)
  for (D in c(1L,3L)) {
    N = 5L; M = 13L
    Y = matrix(rpois(N*M, 1), N, M); Y[1,]=0; Y[,M]=0
    lambda = matrix(rnorm(N*D),N,D)
    dictionary = matrix(rnorm(D*M, sd=15),D,M)
    bound = numeric(N)
    for (i in seq_len(N)) for (m in seq_len(M)) {
      score = lambda[i,] + dictionary[,m]
      log_xi = score - (max(score) + log(sum(exp(score-max(score)))))
      bound[i] = bound[i] + Y[i,m]*sum(exp(log_xi)*(score-log_xi))
    }
    inputs = list(Y)
    if (requireNamespace("Matrix", quietly=TRUE)) inputs[[2]] = Matrix::Matrix(Y,sparse=TRUE)
    for (input in inputs) for (width in c(1L,4L,M)) {
      prepared = .miso_prepare_pass(input,width)
      for (allocated in c(FALSE,TRUE)) {
        result = prepared$pass(prepared$blocks,lambda,dictionary,
                                evaluate_bound=TRUE,compute_allocated=allocated)
        close(result$bound,bound)
        stopifnot(is.null(result$allocated) == !allocated)
      }
    }
  }
  ## Deliberately large eps: LSE evaluates the entropy identity without clipping
  ## log(xi). Floors elsewhere remain unchanged. This is an intentional change.
  prepared = .miso_prepare_pass(matrix(5,1,1),1)
  result = prepared$pass(prepared$blocks,matrix(c(0,-2),1,2),matrix(0,2,1),
                          evaluate_bound=TRUE,compute_allocated=FALSE,eps=0.4)
  close(result$bound,5*log1p(exp(-2)))

  for (S in c(1L,3L)) for (D in c(1L,2L)) for (K in c(1L,4L)) {
    N=5L; M=13L
    Y=matrix(rpois(N*M,1),N,M); Y[1,]=0; Y[,M]=0
    F=.miso_normalize_rows(matrix(rexp(K*M),K,M))
    gamma=.miso_normalize_gamma(array(rexp(S*D*K),c(S,D,K)))
    alpha=array(runif(N*S*D,1,5),c(N,S,D)); beta=matrix(runif(S*D,1,3),S,D)
    omega=.miso_normalize_rows(matrix(rexp(N*S),N,S))
    if (S>1) omega[,S]=0  # An unused motif must have zero scores.
    reference=array(0,c(S,D,K))
    for (s in seq_len(S)) for (i in seq_len(N)) for (m in seq_len(M)) {
      log_weight = vapply(seq_len(D),function(d) digamma(alpha[i,s,d])-log(beta[s,d])+
                            sum(gamma[s,d,]*log(F[,m])),0.0)
      xi = exp(log_weight-max(log_weight)); xi=xi/sum(xi)
      for (d in seq_len(D)) reference[s,d,] = reference[s,d,] +
        omega[i,s]*Y[i,m]*xi[d]*log(F[,m])
    }
    inputs=list(Y,matrix(0,N,M))
    if (requireNamespace("Matrix",quietly=TRUE)) inputs=c(inputs,lapply(inputs,Matrix::Matrix,sparse=TRUE))
    for (input in inputs) for (width in c(1L,4L,M)) {
      expected = if (sum(input)==0) array(0,c(S,D,K)) else reference
      prepared = .miso_prepare_pass(input, width)
      # Deliberately distinct beta0: allocations must use the stored beta.
      prior = matrix(1, S, D)
      streamed = .miso_susie_step(prepared, F, gamma, alpha, beta,
                                  prior, prior, omega, update_F = FALSE)
      dense_counts = .miso_susie_step(prepared, F, gamma, alpha, beta,
                                      prior, prior, omega)
      close(streamed$scores, expected)
      stopifnot(is.null(streamed$C), is.null(dense_counts$scores))
      close(streamed$alpha, dense_counts$alpha)
      close(.miso_gamma_from_scores(streamed$scores, gamma),
            .miso_update_gamma(F, dense_counts$C, gamma))
    }
  }
  if (requireNamespace("Matrix",quietly=TRUE) && isTRUE(capabilities("profmem"))) {
    N=6L; M=10000L; S=11L; D=3L; K=2L
    Y=Matrix::sparseMatrix(i=1:6,j=1:6,x=1,dims=c(N,M))
    F=matrix(1/M,K,M); gamma=array(1/K,c(S,D,K))
    alpha=array(2,c(N,S,D)); beta=matrix(2,S,D); omega=matrix(1/S,N,S)
    prepared=.miso_prepare_pass(Y,100)
    path=tempfile(); on.exit({Rprofmem(NULL);unlink(path)},add=TRUE)
    Rprofmem(path)
    prior=matrix(1,S,D)
    result=.miso_susie_step(prepared,F,gamma,alpha,beta,prior,prior,omega,update_F=FALSE)
    Rprofmem(NULL)
    bytes=suppressWarnings(as.numeric(sub(" .*","",readLines(path))))
    stopifnot(max(bytes,na.rm=TRUE)<S*D*M*8/2,
              identical(dim(result$scores),c(S,D,K)), is.null(result$C))
  }
})
cat("Log-sum-exp and fixed-F streaming tests passed.\n")
