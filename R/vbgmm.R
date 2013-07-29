# Variational Bayesian Gaussian mixture model (VB-GMM)
# X: N x D data matrix
# init: k (1 x 1) or label (1 x n, 1<=label(i)<=k) or center (d x k)
# Reference: Pattern Recognition and Machine Learning by Christopher M. Bishop (P.474)
# Part of the implementation is based on the Matlab code vbgm.m from Michael Chen
# Matlab code: http://www.mathworks.com/matlabcentral/fileexchange/35362-variational-bayesian-inference-for-gaussian-mixture-model
# Author: Yue Li (yueli@cs.toronto.edu)

# Main function vbgmm running four separate functions, namely,
# intialization: initialization of responsibility
# vexp: Variational-Expectation
# vmax: V-Maximimization
# vbound: V-(lower)-bound evaluation
require(pracma)
require(Matrix)


############ bsxfun {pracma} with single expansion (real Matlab style) ############
# expandByRow applies only when x is a square matrix
bsxfun.se <- function(func, x, y, expandByRow=TRUE) {
	
	if(length(y) == 1) return(arrayfun(func, x, y)) else	
		stopifnot(nrow(x) == length(y) || ncol(x) == length(y))
	
	expandCol <- nrow(x) == length(y)	
	expandRow <- ncol(x) == length(y)			
	
	if(expandCol & expandRow & expandByRow) expandCol <- FALSE
	if(expandCol & expandRow & !expandByRow) expandRow <- FALSE
			
	# repeat row (if dim2expand = 1, then length(y) = ncol(x))
	if(expandRow) y.repmat <- matrix(rep(as.numeric(y), each=nrow(x)), nrow=nrow(x))
	
	# repeat col (if dim2expand = 2, then length(y) = nrow(x))
	if(expandCol) y.repmat <- matrix(rep(as.numeric(y), ncol(x)), ncol=ncol(x))
	
	bsxfun(func, x, y.repmat)	
}


############ Compute log(sum(exp(x), dim)) ############
# Compute log(sum(exp(x),dim)) while avoiding numerical underflow
logsumexp <- function(x, margin=1) {
	
	stopifnot(is.matrix(x))
	
	# subtract the largest in each column
	y <- apply(x, margin, max)
					
	x <- bsxfun.se("-", x, y)
	
	s <- y + log(apply(exp(x), margin, sum))
	
	i <- which(!is.finite(s))
	
	if(length(i) > 0) s[i] <- y[i]
	
	s
}


############ Logarithmic Multivariate Gamma function ############
# Compute logarithm multivariate Gamma function.
# Gamma_p(x) = pi^(p(p-1)/4) prod_(j=1)^p Gamma(x+(1-j)/2)
# log Gamma_p(x) = p(p-1)/4 log pi + sum_(j=1)^p log Gamma(x+(1-j)/2)
logmvgamma <- function(x, d) {
	
	s <- size(x)
	
	x <- matrix(as.numeric(x), nrow=1)
			
	x <- bsxfun.se("+", kronecker(matrix(1,d,1), x), (1 - matrix(1:d))/2)
	
	y <- d*(d-1)/4*log(pi) + colSums(lgamma(x))
	
	y <- matrix(as.numeric(y), nrow=s[1], ncol=s[2])
	
	y
}

############ Sort the model  ############
# Sort model paramters in increasing order of averaged means
# of d variables
sort_components <- function(model) {
	
	idx <- order(apply(model$m, 2, mean))
	
	model$m <- model$m[, idx]
	
	model$R <- model$R[, idx]
	
	model$logR <- model$logR[, idx]
		
	model$alpha <- model$alpha[idx]
	
	model$kappa <- model$kappa[idx]
	
	model$v <- model$v[idx]
	
	model$M <- model$M[,,idx]
	
	model
}


############ Vector dot product ############
# handle single row matrix by multiplying each value
# but not sum them up
dot.ext <- function(x,y,mydim) {
	
	if(missing(mydim)) dot(x,y) else {
				
		if(1 %in% size(x) & mydim == 1) x * y else dot(x,y)		
	}
}



############ Main function of VB-GMM ############
# mirprior: miRNA-specific prior

vbgmm <- function(data, init=2, prior, tol=1e-20, maxiter=2e3, mirprior=TRUE, 
	expectedTargetFreq=0.01, verbose=FALSE) {
	
	data <- as.matrix(data)
	
	n <- nrow(data)
	d <- ncol(data)
	
	X <- t(data) # Work with D by N for convenience		
	
	message(sprintf("Running VB-GMM on a %d-by-%d data ...\n", n, d))
	
	if(missing(prior)) {
		
		# miRNA-specific prior with alpha for target-component
		# with defined target frequency (expected much lower
		# than non-target freq)
		if(mirprior & length(init) == 1) {
			
			k <- init
			
			stopifnot(expectedTargetFreq > 0 & expectedTargetFreq < 1)
									
			expectTargetCnt <- round(expectedTargetFreq * n)
			
			backgroundCnt <- rep((n - expectTargetCnt)/(k-1), k-1)
			
			prior <- list(
				alpha = c(expectTargetCnt, backgroundCnt),
				kappa = 1,
				m = as.matrix(rowMeans(X)),
				v = d+1,
				M = diag(1,d,d) # M = inv(W)
			)
			
		} else { # more general prior with equal alpha
				
			prior <- list(
				alpha = rep(1,k),
				kappa = 1,
				m = as.matrix(rowMeans(X)),
				v = d+1,
				M = diag(1,d,d) # M = inv(W)
			)
		}
	} else {
		
		stopifnot(
			all(names(prior) %in% c("alpha","kappa","m","v","M")) &
			all(sapply(prior, is.numeric)) & nrow(prior$m) == d &
			ncol(prior$m) == 1 &
			nrow(prior$M) == d & ncol(prior$M) == d)					
	}
	
	if(mirprior & length(init) != 1) {
		warning("mirprior is TRUE but init is not scalar (k)! Proceed as if mirprior were FALSE") }
	
	# lower variational bound (objective function)
	L <- rep(-Inf, maxiter)	
	converged <- FALSE	
	t <- 1
	
	model <- list()
	
	model$R <- initialization(X, init) # initialize responsibility R
	
	while(!converged & t < maxiter) {				
		t <- t + 1
		model <- vmax(X, model, prior)
		model <- vexp(X, model)
		L[t] <- vbound(X, model, prior)/n
		converged <- abs(L[t] - L[t-1]) < tol * abs(L[t])
		
		if(verbose) message(sprintf("VB-EM-%d: L = %.6f", t, L[t]))
	}
	
	L <- L[2:t]
	
	model <- sort_components(model)
	
	label <- rep(0, n)
	
	label <- apply(model$R, 1, which.max)		
	
	# unique to have consecutive label eg 2,3,6 changed to 1,2,3
	# label <- match(label, sort(unique(label)))	
	
	if(converged) message(sprintf("Converged in %d steps.\n", t-1)) else
	warnings(sprintf("Not converged in %d steps.\n", maxiter))
	
	list(label=label, R=model$R, mu=model$m, full.model=model, L=L)
}


############ Initialization of responsibility (intialization) ############
initialization <- function(X, init) {
	
	d <- nrow(X)
		
	n <- ncol(X)
	
	stopifnot(length(init) %in% c(1, n) ||
		(nrow(init) == d  & ncol(init) == k))
		
	if(length(init) == 1) { # init = k gaussian components
		
		k <- init
		
		idx <- sample(1:n, k)
		
		m <- X[,idx,drop=F]		
		
		label <- apply(bsxfun.se("-", crossprod(m, X),
			as.matrix(dot.ext(m,m,1)/2)), 2, which.max)
		
		# unique to have consecutive label eg 2,3,6 changed to 1,2,3
		label <- match(label, sort(unique(label)))
		
		while(k != length(unique(label))) {
			
			idx <- sample(1:n, k)
		
			m <- X[,idx,drop=F]
			
			label <- apply(bsxfun.se("-", crossprod(m, X),
				as.matrix(dot.ext(m,m,1)/2)), 2, which.max)
				
			label <- match(label, unique(label))			
		}
		
		R <- as.matrix(sparseMatrix(1:n, label, x=1))
		
	} else {
	
		if(length(init) == n) { # initialize with labels
			
			label <- init		
			k <- max(label)		
			R <- as.matrix(sparseMatrix(1:n, label, x=1))
			
		} else {
		
			if(!is.null(dim(init))) {
				
				if(nrow(init) == d  & ncol(init) == k) { # initialize with centers
				
					k <- ncol(init)		
					m <- init		
					label <- apply(bsxfun.se("-", crossprod(m, X),
						as.matrix(dot.ext(m,m,1)/2)), 2, which.max)		
					R <- as.matrix(sparseMatrix(1:n, label, x=1))
					
				} else stop(message("Invalid init."))				
			}
		}
	}

	R	
}


############ Variational-Maximimization ############
vmax <- function(X, model, prior) {
					
	alpha0 <- prior$alpha	
	kappa0 <- prior$kappa
	m0 <- prior$m;
	v0 <- prior$v;
	M0 <- prior$M;
	R <- model$R;

	nk <- apply(R, 2, sum) # 10.51
	alpha <- alpha0 + nk # 10.58
	nxbar <- X %*% R
	kappa <- kappa0 + nk # 10.60		
	m <- bsxfun.se("*", bsxfun.se("+", nxbar, kappa0 * m0), 1/kappa) # 10.61
	v <- v0 + nk # 10.63 (NB: no 1 in the matlab code)
	
	d <- nrow(m)
	k <- ncol(m)
	
	M <- array(0, c(d, d, k))
	sqrtR <- sqrt(R)
	
	xbar <- bsxfun.se("*", nxbar, 1/nk) # 10.52
	xbarm0 <- bsxfun.se("-", xbar, m0)
	w <- (kappa0 * nk) * (1/(kappa0 + nk))
	
	for(i in 1:k) {
		
		Xs <- bsxfun.se("*", bsxfun.se("-", X, xbar[,i]), t(sqrtR[,i]))
		
		xbarm0i <- xbarm0[,i]
		
		# 10.62
		M[,,i] <- M0 + tcrossprod(Xs, Xs) + w[i] * tcrossprod(xbarm0i, xbarm0i)
	}		
	
	model$alpha <- alpha
	model$kappa <- kappa
	model$m <- m
	model$v <- v
	model$M <- M # Whishart: M = inv(W)
	
	model
}


############ Variational-Expectation ############
vexp <- function(X, model) {
	
	alpha <- model$alpha	# Dirichlet
	kappa <- model$kappa	# Gaussian
	m <- model$m			# Gasusian
	v <- model$v			# Whishart
	M <- model$M			# Whishart: inv(W) = V'*V
	
	n <- ncol(X)
	d <- nrow(m)
	k <- ncol(m)
	
	logW <- array(0, dim=c(1,k))
	EQ <- array(0, dim=c(n,k))
	
	for(i in 1:k) {
		
		U <- chol(M[,,i])
		
		logW[i] <- -2 * sum(log(diag(U)))
		
		Q <- solve(t(U), bsxfun.se("-", X, m[,i]))
		
		EQ[,i] <- d/kappa[i] + v[i] * dot.ext(Q,Q,1)	# 10.64
	}
	
	vd <- bsxfun.se("-", matrix(rep(v+1, d),nrow=d,byrow=T), as.matrix(1:d))/2
	
	ElogLambda <- colSums(digamma(vd)) + d*log(2) + logW	# 10.65	
	Elogpi <- digamma(alpha) - digamma(sum(alpha))			# 10.66
	
	logRho <- (bsxfun.se("-", EQ, 2*Elogpi + ElogLambda - d*log(2*pi)))/(-2) # 10.46
	logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
	R <- exp(logR)
	
	model$logR <- logR
	model$R <- R
	
	model
}



############ Variational-(lower)-Bound Evaluation ############
vbound <- function(X, model, prior) {
	
	alpha0 <- prior$alpha
	kappa0 <- prior$kappa
	m0 <- prior$m
	v0 <- prior$v
	M0 <- prior$M

	alpha <- model$alpha	# Dirichlet
	kappa <- model$kappa	# Gaussian
	m <- model$m			# Gasusian
	v <- model$v			# Whishart
	M <- model$M			# Whishart: inv(W) = V'*V
	R <- model$R
	logR <- model$logR
	
	d <- nrow(m)
	k <- ncol(m)
	nk <- colSums(R)									# 10.51
	
	Elogpi <- digamma(alpha) - digamma(sum(alpha))		# 10.66
	
	Epz = dot(nk, Elogpi)								# 10.72
	Eqz = dot(as.numeric(R), as.numeric(logR))			# 10.75
	# logCalpha0 = lgamma(k * alpha0) - k * lgamma(alpha0) # for scalar alpha0
	logCalpha0 = lgamma(sum(alpha0)) - sum(lgamma(alpha0))
	# Eppi <- logCalpha0+(alpha0-1)*sum(Elogpi) # for scalar alpha0
	Eppi <- logCalpha0+dot(alpha0-1, Elogpi)			# 10.73		
	logCalpha <- lgamma(sum(alpha))-sum(lgamma(alpha))
	Eqpi = dot(alpha-1, Elogpi) + logCalpha				# 10.76
	
	# part of 10.70
	L <- Epz - Eqz + Eppi - Eqpi
	
	U0 <- chol(M0)
	sqrtR <- sqrt(R)
	xbar <- bsxfun.se("*", X %*% R, 1/nk)				# 10.52
	
	logW <- array(0, dim = c(1, k))
	trSW <- array(0, dim = c(1, k))
	trM0W <- array(0, dim = c(1, k))
	xbarmWxbarm <- array(0, dim = c(1, k))
	mm0Wmm0 <- array(0, dim = c(1, k))
	
	for(i in 1:k) {
		
		U <- chol(M[,,i])
		logW[i] <- -2 * sum(log(diag(U)))
		
		Xs <- bsxfun.se("*", bsxfun.se("-", X, as.matrix(xbar[,i,drop=F])), t(sqrtR[,i,drop=F]))
		V <- chol(tcrossprod(Xs, Xs)/nk[i])
		Q <- solve(U, V)
		# equivalent to tr(SW)=trace(S/M)
		trSW[i] <- dot(as.numeric(Q), as.numeric(Q))
		Q <- solve(U, U0)
		trM0W[i] <- dot(as.numeric(Q), as.numeric(Q))
		
		q <- solve(t(U), xbar[,i,drop=F]-m[,i,drop=F])
		xbarmWxbarm[i] = dot(q, q)
		q <- solve(t(U), m[,i,drop=F]-m0)
		mm0Wmm0[i] <- dot(q, q)
	}
	
	vd <- bsxfun.se("-", matrix(rep(v+1, d),nrow=d,byrow=T), as.matrix(1:d))/2	
	ElogLambda <- colSums(digamma(vd)) + d*log(2) + logW	# 10.65
	
	# first half of 10.74
	Epmu <- sum(d*log(kappa0/(2*pi))+ElogLambda-d*kappa0/kappa-kappa0*(v*mm0Wmm0))/2
	logB0 <- v0*sum(log(diag(U0)))-0.5*v0*d*log(2)-logmvgamma(0.5*v0,d) # B.79
	# second half of 10.74	
	EpLambda <- k*logB0+0.5*(v0-d-1)*sum(ElogLambda)-0.5*dot(v,trM0W)
	
	
	Eqmu <- 0.5*sum(ElogLambda+d*log(kappa/(2*pi)))-0.5*d*k	# 10.77 (1/2)
	logB <- (-v) * (logW+d*log(2))/2 - logmvgamma(0.5*v, d)	# B.79
	EqLambda <- 0.5*sum((v-d-1)*ElogLambda-v*d)+sum(logB) 	# 10.77 (2/2)
	
	EpX <- 0.5*dot(nk, ElogLambda-d/kappa-v*trSW-v*xbarmWxbarm-d*log(2*pi))	# 10.71
		
	L <- L+Epmu-Eqmu+EpLambda-EqLambda+EpX	# 10.70
	
	L
}








