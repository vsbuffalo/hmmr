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
        hstates[i] <- sample(states, 1,  prob=hmm@A[hstates[i-1], ])
        seq[i] <- sample(emissions, 1, prob=hmm@B[hstates[i-1], ])
    }

    if (!is.null(dimnames(hmm@B))) {
        # make results labelled factor if emissoin probability matrix
        # has labels
        seq <- factor(seq, labels=colnames(hmm@B))
        hstates <- factor(hstates, labels=rownames(hmm@B))
    }
    list(sequence=seq, states=hstates)
})

setMethod("forward", signature=c(hmm="HMM", x="numeric"), function(hmm, x) {
    # Forward algorithm, with scaling for numerical stability.
    states <- seq_len(nrow(hmm@A)) # number of states
    L <- length(x)
    alpha <- matrix(0, nrow=nrow(hmm@A), ncol=L)
    scaling <- numeric(L)

    # main recursion component
    for (i in seq_along(x)) {
        for (k in states) {
            if (i == 1) {
                # initial state is simply initial probs * emissions
                # for a state and observed symbol
                alpha_sum <- hmm@A0[k] * hmm@B[k, x[i]]
            } else {
                # marginalize over product of last state's
                # probabilities and transition to this state
                alpha_sum <- sum(last_alpha * hmm@A[, k])
            }
            alpha[k, i] <- alpha_sum * hmm@B[k, x[i]]
        }
        last_alpha <- alpha[, i]
    }
    return(alpha)
})

hmm <- HMM(swaps_trans, coin_probs, initial_probs)
flips <- simulate(hmm, 10)

apply(forward(hmm, flips$sequence), 2, which.max)
apply(backward(hmm, flips$sequence), 2, which.max)
flips$sequence
flips$state


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
    return(beta)
})


forwardBackward <-
function(sequence, transitions, emission_probs, initial_probs) {
    fwd <- forward(sequence, transitions, emission_probs, initial_probs)
    bwd <- backward(sequence, transitions, emission_probs, initial_probs)

    p_fwd <- sum(transitions[, sequence[length(sequence)]], fwd[, ncol(fwd)])
    p_bwd <- sum(initial_probs * emission_probs[, sequence[1]] * bwd[, 1])

    posterior <- fwd * bwd / p_fwd
    return(posterior)
}



### Coin Flip Simulation
## simulated biased coin
swaps_trans <- matrix(c(0.98, 0.02, 0.02, 0.98), nrow=2, ncol=2)
# H, T
initial_probs <- c(0.1, 0.9)
# (fair, biased) x (p(H), p(T))
coin_probs <- matrix(c(0.5, 0.8, 0.5, 0.2), nrow=2, ncol=2,
                     dimnames=list(c("F", "B"), c("H", "T")))

hmm <- HMM(swaps_trans, coin_probs, initial_probs)
flips <- simulate(hmm, 10)
