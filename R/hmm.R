## hmm.R -- Hidden Markov Models

# General notes: states are held in indices of transition matrices and
# emission index

setClass("HMM",
         representation(A="matrix",
                        B="matrix",
                        A0="numeric"))

HMM <-
# Constructor method    
function(A, B, A0, states=NULL, symbols=NULL) {
    stopifnot(nrow(A) == ncol(A))
    stopifnot(all(colSums(A) == 1))
    stopifnot(all(rowSums(B) == 1))

    # label matrices 
    if (!is.null(states))
        rownames(A) <- colnames(A) <- rownames(B) <- names(A0) <- states
    if (!is.null(symbols))
        colnames(B) <- symbols
    new("HMM", A=A, B=B, A0=A0)
}

setGeneric("forward", function(hmm, x) standardGeneric("forward"))
setGeneric("backward", function(hmm, x) standardGeneric("backward"))
setGeneric("forwardBackward", function(hmm, x) standardGeneric("forwardBackward"))

setMethod("show", "HMM", function(object) {
    cat(sprintf("HMM Object with %d symbols and %d states\n", nrow(object@B), ncol(object@A)))
    cat("transitions:\n")
    print(object@A)
    cat("emissions:\n")
    print(object@B)
    cat("initial probabilities:\n")
    print(object@A0)
})

setMethod("simulate", signature=c(object="HMM"),
    function(object, nsim, seed=NULL, ...) {
    # Function that simulates observed sequences from an HMM model of a
    # transition matrix, initial state probabilities, and emission
    # probabilities.
    set.seed(seed)
    states <- seq_len(nrow(hmm@A))
    emissions <- seq_len(ncol(hmm@B))
        
    # generate n states from the transition matrix
    hstates <- integer(nsim)
    seq <- integer(nsim)
    hstates[1] <- sample(states, 1, prob=hmm@A0)
    seq[1] <- sample(emissions, 1, hmm@A0[hstates[1]])
    for (i in seq(2, nsim)) {
        # sample a state, using the last state's probabilities
        hstates[i] <- sample(states, 1, prob=hmm@A[hstates[i-1], ])
        seq[i] <- sample(emissions, 1, prob=hmm@B[hstates[i-1], ])
    }

    if (!is.null(dimnames(hmm@B))) {
        # make results labelled factor if emissoin probability matrix
        # has labels
        seq <- factor(seq, levels=states, labels=colnames(hmm@B))
        hstates <- factor(hstates, levels=states, labels=rownames(hmm@B))
    }
    list(sequence=seq, states=hstates)
})

setMethod("forward", signature=c(hmm="HMM", x="numeric"), function(hmm, x) {
    # Forward algorithm, with scaling for numerical stability.
    states <- seq_len(nrow(hmm@A)) # number of states
    L <- length(x)
    alpha <- matrix(0, nrow=nrow(hmm@A), ncol=L)
    scaling <- numeric(L)
    last_alpha <- alpha[, 1] <- hmm@A0# * hmm@B[, x[1]]

    # main recursion component
    for (i in seq_along(x)) {
        for (k in states) {
            alpha[k, i] <- sum(last_alpha * hmm@A[, k]) * hmm@B[k, x[i]]
        }
        last_alpha <- alpha[, i]
    }
    dimnames(alpha) <- list(rownames(hmm@B), seq_along(x))
    return(alpha)
})

setMethod("backward", signature=c(hmm="HMM", x="numeric"), function(hmm, x) {
    # Backward algorithm, with scaling for numerical stability.
    states <- seq_len(nrow(hmm@A))
    L <- length(x)
    beta <- matrix(0, nrow=nrow(hmm@A), ncol=L)
    last_beta <- beta[, L] <- rep(1, length(states))

    # main recursion component
    for (i in seq(L-1, 1)) {
        for (k in states) {
            beta[k, i] <- sum(hmm@B[, x[i+1]] * last_beta * hmm@A[k, ])
        }
        last_beta <- beta[, i]
    }
    dimnames(beta) <- list(rownames(hmm@B), seq_along(x))
    return(beta)
})

setMethod("forwardBackward", signature=c("HMM", x="numeric"), function(hmm, x) {
    fwd <- forward(hmm, x)
    bwd <- backward(hmm, x)
    p_fwd <- fwd[, length(x)]
    posterior <- (fwd * bwd) / sum(p_fwd)
    return(posterior)
})
