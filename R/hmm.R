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
    stopifnot(all(colSums(A)) == 1)
    stopifnot(all(rowSums(B)) == 1)

    # label matrices 
    if (!is.null(states))
        rownames(A) <- colnames(A) <- rownames(B) <- names(A0) <- states
    if (!is.null(symbols))
        colnames(B) <- symbols
    new("HMM", A=A, B=B, A0=A0)
}

setGeneric("forward")
setGeneric("backward")
setGeneric("forwardBackward")

setMethod("show", "HMM", function(x) {
    cat(sprintf("HMM Object with %d symbols and %d states\n", nrow(x@B), ncol(x@A)))
    cat("transitions:\n")
    print(x@A)
    cat("emissions:\n")
    print(x@B)
    cat("initial probabilities:\n")
    print(x@A0)
})

forward <-
# Forward algorithm for HMMs - O(nm^2)
#
# Computes the forward probability, the probability of emitting an
# observed sequence and reaching a particular state.
function(sequence, transitions, emission_probs, initial_probs) {
    states <- seq_len(nrow(transitions))

    # m X n where m is number of states and n is the number of
    # observations.
    m <- length(states)
    n <- length(sequence)
    
    # out forward probabilities over each observation
    fwd_probs <- matrix(0, nrow=m, ncol=n)
    
    last_fwd_probs <- initial_probs

    # main recursion component
    for (i in seq_along(sequence)) {
        x_i <- sequence[i]
        for (l in states) {
            if (i == 1) {
                fwd_prob_sum <- last_fwd_probs[l]
            } else {
                fwd_prob_sum <- sum(last_fwd_probs * transitions[, l])
            }
            fwd_probs[l, i] <- fwd_prob_sum * emission_probs[l, x_i]
        }
        last_fwd_probs <- fwd_probs[, i]
    }
    return(fwd_probs)
}

## x <- forward(flips$sequence, swaps_trans, coin_probs, initial_probs)
## apply(x, 2, which.max)

backward <-
# Backward algorithm for HMMs
#
# Computes the backward probability, the probability of being a state
# at position i, and emitting the observed sequenced x_i, x_{i+1},
# ... x_n.
function(sequence, transitions, emission_probs, initial_probs) {
    states <- seq_len(nrow(transitions))

    # m X n where m is number of states and n is the number of
    # observations.
    m <- length(states)
    n <- length(sequence)
    
    # out forward probabilities over each observation
    bwd_probs <- matrix(0, nrow=m, ncol=n)
    # main recursion component
    for (i in seq(n, 1)) {
        for (l in states) {
            if (i == n) {
                bwd_probs[, i] <- rep(1, length(states))    
            } else {
                tmp <- sum(emission_probs[, sequence[i+1]]*last_bwd_probs*transitions[l, ])
                bwd_probs[l, i] <- tmp
            }
        }
        last_bwd_probs <- bwd_probs[, i]
    }
    return(bwd_probs)
}
    
forwardBackward <-
function(sequence, transitions, emission_probs, initial_probs) {
    fwd <- forward(sequence, transitions, emission_probs, initial_probs)
    bwd <- backward(sequence, transitions, emission_probs, initial_probs)

    p_fwd <- sum(transitions[, sequence[length(sequence)]], fwd[, ncol(fwd)])
    p_bwd <- sum(initial_probs * emission_probs[, sequence[1]] * bwd[, 1])

    posterior <- fwd * bwd / p_fwd
    return(posterior)
}
