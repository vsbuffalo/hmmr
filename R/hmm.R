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
setGeneric("backward", function(hmm, x, scaling) standardGeneric("backward"))
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
    last_alpha <- alpha[, 1] <- hmm@A0

    # main recursion component
    for (i in seq_along(x)) {
        for (k in states) {
            alpha[k, i] <- sum(last_alpha * hmm@A[, k]) * hmm@B[k, x[i]]
        }
        scaling[i] <- c_i <- sum(alpha[, i])
        last_alpha <- alpha[, i] <- alpha[, i]/c_i
    }
    dimnames(alpha) <- list(rownames(hmm@B), seq_along(x))
    return(list(alpha=alpha, scaling=1/scaling))
})

setMethod("backward", signature=c(hmm="HMM", x="numeric"), function(hmm, x, scaling) {
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
        beta[, i] <- last_beta <- scaling[i]* beta[, i]
    }
    dimnames(beta) <- list(rownames(hmm@B), seq_along(x))
    return(beta)
})

setMethod("forwardBackward", signature=c("HMM", x="numeric"), function(hmm, x) {
    fwd <- forward(hmm, x)
    bwd <- backward(hmm, x, fwd$scaling)
    alpha_beta <- fwd$alpha * bwd
    posterior <- sweep(alpha_beta, 2, colSums(alpha_beta), "/")
    return(posterior)
})

plotPosterior <- function(states, posterior) {
    postdf <- data.frame(t(posterior))
    postdf$max <- factor(apply(postdf, 1, which.max), labels=rownames(posterior))
    postdf$i <- seq_len(ncol(posterior))
    postdf$state <- states
    p <- ggplot(postdf) + geom_line(aes_string(x="i", y=colnames(postdf)[1])) + geom_rug(aes(x=i, color=state))
    p <- p + ylab(sprintf("posterior probability of %s", colnames(postdf)[1]))
    p + xlab("index")
}
