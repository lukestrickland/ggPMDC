### Hierarchical tools ------------------------------------------------
##' Get a n-cell matrix
##'
##' Constructs a matrix, showing how many responses to in each
##' cell. The function checks whether the format of \code{n} and \code{ns}
##' conform.
##'
##' \code{n} can be:
##' \enumerate{
##' \item an integer for a balanced design,
##' \item a matrix for an unbalanced design, where rows are subjects and
##' columns are cells. If the matrix is a row vector, all subjects
##' have the same \code{n} in each cell. If it is a column vector, all
##' cells have the same \code{n}. Otherwise each entry specifies the \code{n}
##' for a particular subject x cell combination. See below for concrete
##' examples.}
##'
##' @param model a model object.
##' @param n number of trials.
##' @param ns number of subjects.
##' @examples
##' model <- BuildModel(
##'   p.map     = list(A = "1", B = "R", t0 = "1", mean_v = "M", sd_v = "M",
##'                   st0 = "1"),
##'   match.map = list(M = list(s1 = 1, s2 = 2)),
##'   constants = c(sd_v.false = 1, st0 = 0),
##'   factors   = list(S = c("s1","s2")),
##'   responses = c("r1", "r2"),
##'   type      = "norm")
##'
##' #######################30
##' ## Example 1
##' #######################30
##' GetNsim(model, ns = 2, n = 1)
##' #      [,1] [,2]
##' # [1,]    1    1
##' # [2,]    1    1
##'
##' #######################30
##' ## Example 2
##' #######################30
##' n <- matrix(c(1:2), ncol = 1)
##' #      [,1]
##' # [1,]    1  ## subject 1 has 1 response for each cell
##' # [2,]    2  ## subject 2 has 2 responses for each cell
##'
##' GetNsim(model, ns = 2, n = n)
##' #      [,1] [,2]
##' # [1,]    1    1
##' # [2,]    2    2
##'
##' #######################30
##' ## Example 3
##' #######################30
##' n <- matrix(c(1:2), nrow = 1)
##' #      [,1] [,2]
##' # [1,]    1    2
##' GetNsim(model, ns = 2, n = n)
##' #     [,1] [,2]
##' # [1,]   1    2 ## subject 1 has 1 response for cell 1 and 2 responses for cell 2
##' # [2,]   1    2 ## subject 2 has 1 response for cell 1 and 2 responses for cell 2
##'
##' #######################30
##' ## Example 4
##' #######################30
##' n <- matrix(c(1:4), nrow=2)
##' #      [,1] [,2]
##' # [1,]    1    3
##' # [2,]    2    4
##' ggdmc::GetNsim(model, ns = 2, n = n)
##' #      [,1] [,2]
##' # [1,]    1    3 ## subject 1 has 1 response for cell 1 and 3 responses for cell 2
##' # [2,]    2    4 ## subject 2 has 2 responses for cell 1 and 4 responses for cell 2
##
##' @export
GetNsim <- function (model, n, ns)
## Use only in simulate_many R function
{
  # ns <- 2
  # n <- 1
  if (ns <= 1) stop("Use simulate.model instead to simulate one participant.")
  if (is.vector(n) & (length(n) != 1))
  {
    if (!is.matrix(n)) stop("n must be a scalar, a vector, or a matrix")
  }

  facs  <- attr(model, "factors")
  ## sapply(facs, length) return the numbers of level in each factor in a named
  ## vector. Then prod multiplies them to get the number of cell
  ncell <- prod(sapply(facs, length))

  ## Untested
  if (is.matrix(n))
  {
    dim1 <- dim(n)[1]
    dim2 <- dim(n)[2]
    if (dim1 == 1) {
      nmat <- matrix(rep(n, each = ns), ns) ## only cells differ
    } else if (dim2 == 1) {
      nmat <- matrix(rep.int(n, ncell), ns) ## only subjects differ
    } else {
      nmat <- n
    }
  } else {
    nmat <- matrix(rep.int(n, ns*ncell), ns)
  }

  dim1 <- dim(nmat)[1]   ## n has been altered in ifelse
  dim2 <- dim(nmat)[2]

  if ( ns != dim1 )    stop(paste0("The n matrix must have ", ns, " rows"))
  if ( ncell != dim2 ) stop(paste0("The n matrix must have ", ncell, " columns"))
  return(nmat)
}

ismanymodels <- function(model, ns = NA)
## Used in simulate_many
{
  if (is.list(model))
  {
    if (length(model) != ns)
      stop("number of participants not equal to number of models")
    out <- TRUE
  }
  else
  {
    if (is.na(ns)) stop("Must indicate the number of participants")
    out <- FALSE
  }
  return(out)
}

##' Extract trial log likelihoods
##'
##' This function simply run trial_loglik to loop through one subject after
##' another to extracts trial_log_likes from a list of subject fits and
##' concatanates the result into an array.
##'
##' @param samples posterior samples
##' @param thin thinnng length
##' @param verbose whether print information
##' @export
trial_loglik_hier <- function(samples, thin = 1, verbose=FALSE)
{
  check <- function(x) {
    nmc <- min(x[1,])

    if ( !all(x[3,1]==x[3,-1]) )
      warning(paste("Subjects do not all have the same number of interations, using first",
                    nmc,"for all."))

    if ( !all(x[2,1]==x[2,-1]) )
      stop("Subjects must have the same number of chains")
  }
  PrintSize <- function(x, thin) {

    tll <- trial_loglik(x[[1]], thin)

    nsub <- length(x)
    sdim <- dim(tll); ## sdim ##  ntrial x nchain  x nnmc_thin
    size <- sum(nsub * prod(sdim))

    ## nchain and nmc per subject
    nmc_nchain_mat  <- sapply(x, function(x) { dim(x$log_likelihoods) } )

    dimnames(nmc_nchain_mat) <- list(c("Chains", "Trials"), names(x))
    cat("Log-likelihood dimension\n")
    print(nmc_nchain_mat)

    names(sdim) <- c("Trials","Chains","Iterations")
    cat("\nSubject 1\n")
    print(sdim)
    cat("Total log-likelihoods:", round(size/1e6,2), "millions)\n")
    return(NULL)
  }

  if(verbose) PrintSize(samples, thin)

  tlls <- lapply(samples, ggdmc:::trial_loglik, thin)

  sdims <- sapply(tlls, dim)  ## DMC nmc_thin, nchain, ntrial
  nmc <- min(sdims[3,])       ## nnmc_thin
  check(sdims)

  ### nmc_thin x nchain x (ntrial x nsub)
  out <- array(dim=c(nmc, dim(tlls[[1]])[2], sum(sdims[1,])))
  nsub <- length(samples)

  start <- 1;
  end <- sdims[1,1]

  for (i in 1:nsub)
  {
    out[,,start:end] <- aperm(tlls[[i]], c(3, 2, 1))
    if (i < nsub)
    {
      start <- end+1
      end   <- start - 1 + sdims[1, i+1]
    }
  }
  return(out)
}

### MCMC ------------------------------------------------
##' Create a MCMC list
##'
##' This funciton is used by ConvertChains function.
##'
##' @param x posterior samples
##' @param start start from which iteration
##' @param end end at which iteration
##' @param pll a Boolean switch for extract posterior log-likelihoods
##' @param ll a Boolean switch to extract log likelihoods
##' @export
mcmc_list.model <- function(x, start = 1, end = NA, pll = TRUE, ll = FALSE)
{

  if (is.null(x$theta)) stop("Use hyper mcmc_list")
  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")
  if (pll) {

    message("Convert to posterior log-likelihood")
    lp <- x$summed_log_prior[,start:end] + x$log_likelihoods[,start:end]
    colnames(lp) <- 1:ncol(lp)

    step1 <- lapply(data.frame(t(lp)), function(xx) {
      coda::mcmc(as.matrix(xx), start, end, thin = 1)
    })
    out <- coda::mcmc.list(step1) ## log-posterior likelihood


  } else if (ll) {
    message("Convert to log-likelihood")
    ll <- x$log_likelihoods[,start:end]
    colnames(ll) <- 1:ncol(ll)

    step1 <- lapply(data.frame(t(ll)), function(xx) {
      coda::mcmc(as.matrix(xx), start, end, thin = 1)
    })
    out <- coda::mcmc.list(step1) ## log-posterior likelihood

  } else {
    message("Convert to theta")
    # lp <- x$log_likelihoods[,start:end]
    out <- theta2mcmclist(x, start = start, end = end, thin = 1)
  }

  attr(out, "nchain") <- x$n.chain
  attr(out, "npar")   <- x$n.pars
  attr(out, "thin")   <- x$thin
  attr(out, "iter")   <- start:end
  attr(out, "pnames") <- x$p.names
  attr(out, "nmc")    <- x$nmc
  attr(out, "start")  <- start
  attr(out, "end")    <- end
  return(out)
}

##' Prepare posterior samples for plotting functions version 1
##'
##' Convert MCMC chains to a data frame for plotting functions. This funciton is
##' used by autocorr.
##'
##' @param x posterior samples
##' @param start which iteration to start
##' @param end end at which iteration
##' @param pll a Boolean switch to extract posterior log likelihoods
##' @param ll a Boolean switch to extract log likelihoods
##' @export
ConvertChains <- function(x, start = 1, end = NA, pll = TRUE, ll = FALSE)
{
  if (x$n.chains == 1) stop ("MCMC needs multiple chains to check convergence")
  if (is.null(x$theta)) stop("Use hyper mcmc_list")
  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")

  nchain <- x$n.chain
  npar   <- x$n.pars
  pnames <- x$p.names
  iter   <- start:end
  mcmclist <- mcmc_list.model(x, start, end, pll, ll)

  v <- lapply(seq_along(mcmclist), function(k) {
    tmp1 <- sapply(mcmclist[k][[1]], c)   ## concatenate function

    if (pll) {
      d <- data.frame(Iteration = iter,Parameter = "lp", value = tmp1)
    } else if (ll) {
      d <- data.frame(Iteration = iter,Parameter = "ll", value = tmp1)
    }
    else {
      d <- data.frame(Iteration = rep(iter, npar),
        Parameter = rep(pnames, each = length(iter)), value = tmp1)
    }

    d$Chain <- k
    d[, c("Iteration", "Chain", "Parameter", "value")]
  })


  D <- data.table::rbindlist(v)
  D$Parameter <- factor(D$Parameter)
  D$Chain <- factor(D$Chain)
  attr(D, "nThin")       <- x$thin
  attr(D, "nIterations") <- end
  attr(D, "nChains")     <- nchain
  attr(D, "nParameters") <- npar

  return(D)
}


ConvertChains2 <- function(x, pll)
{
  nchain <- attr(x, "nchain")
  npar <- attr(x, "npar")
  pnames <- attr(x, "pnames")
  start <- attr(x, "start")
  end <- attr(x, "end")

  if (nchain == 1) stop ("MCMC needs multiple chains to check convergence")
  # if (is.null(x$theta)) stop("Use hyper mcmc_list")
  if ( is.na(end) ) end <- attr(x, "nmc")
  if ( end <= start ) stop("End must be greater than start")
  iter   <- start:end

  v <- lapply(seq_along(x), function(k) {
    tmp1 <- sapply(x[k][[1]], c)
    if (pll) {
      d <- data.frame(Iteration = iter,
        Parameter = "lp", value = tmp1)
    } else {
      d <- data.frame(Iteration = rep(iter, npar),
        Parameter = rep(pnames, each = length(iter)), value = tmp1)
    }

    d$Chain <- k
    d[, c("Iteration", "Chain", "Parameter", "value")]
  })


  D <- data.table::rbindlist(v)
  D$Parameter <- factor(D$Parameter)
  D$Chain <- factor(D$Chain)
  attr(D, "nThin")       <- attr(x, "thin")
  attr(D, "nIterations") <- end
  attr(D, "nChains")     <- nchain
  attr(D, "nParameters") <- npar

  return(D)
}

##' Calculate the autocorrelation of a vector
##'
##' Calculate the autocorrelation of a vector.
##'
##' Internal function used by \code{\link{ggs_autocorrelation}}.
##'
##' @param x a vector storing parameter values
##' @param nLags the maximum number of lags
##' @return A data.frame
##' @export
##' @examples
##' ## Calculate the autocorrelation of a simple vector
ac <- function(x, nLags = 50)
{
  tmp <- ac_(x, nLags)
  return(  list(Lag = 1:nLags, Autocorrelation=tmp[,1]))
}


### Generic  ---------------------------------------------
##' Retrieve information of operating system
##'
##' A wrapper function to extract system information from \code{Sys.info}
##' and \code{.Platform}
##'
##' @examples
##' get_os()
##' ## sysname
##' ## "linux"
##' @export
get_os <- function()
{
  sysinf <- Sys.info()
  ostype <- .Platform$OS.type

  ## Probe using Sys.info: Windows, Linux or Darwin
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") os <- "osx"
  } else {
    ## If something gets wrong with Sys.info, probe using .Platform
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))   os <- "osx"
    if (grepl("unix", R.version$os))      os <- "osx"
    if (grepl("linux-gnu", R.version$os)) os <- "linux"
    if (grepl("mingw32", R.version$os))   os <- "windows"
  }
  tolower(os)
}


# ac <- function (x, nlag) {
#   out <- data.frame(Lag = 1:nlag, Autocorrelation = cor(ac_(x, nlag),
#     use = "pairwise.complete.obs")[, 1])
#   return(out)
# }


##' Extract parameter names from a model object
##'
##' @param x a model object
##'
##' @export
GetPNames <- function(x) { return(names(attr(x, "p.vector"))) }

checklba <- function(x)
{
  model <- attr(x, "model")
  if (attr(model, "type") == "norm" )
  {
    parnames <- attr(model, "par.names")
    if ( (which(parnames == "A")      != 1) |
         (which(parnames == "B")      != 2) |
         (which(parnames == "t0")     != 3) |
         (which(parnames == "mean_v") != 4) |
         (which(parnames == "sd_v")   != 5) |
         (which(parnames == "st0")    != 6) )
    {
      cat("Your p.vector is order as: ", parnames, "\n")
      message("It must be in the order of: A, B, t0, mean_v, sd_v, & st0.")
      stop("Check p.map")
    }
  }

}


