## simulate.R -- simulate sequences from HMM setup

simulateSequence <-
# Function that simulates observed sequences from an HMM model of a
# transition matrix, initial state probabilities, and emission
# probabilities.
function(n, transitions, emission_probs, initial_probs) {
    stopifnot(nrow(transitions) == ncol(transitions))
    stopifnot(nrow(transitions) == nrow(emission_probs))
    stopifnot(nrow(transitions) == length(initial_probs))
   
    states <- seq_len(nrow(transitions))
    emissions <- seq_len(ncol(emission_probs))
        
    # generate n states from the transition matrix
    hstates <- integer(n)
    seq <- integer(n)
    hstates[1] <- sample(states, 1, prob=initial_probs)
    seq[1] <- sample(emissions, 1, emission_probs[hstates[1]])
    for (i in 2:n) {
        # sample a state, using the last state's probabilities
        hstates[i] <- sample(states, 1,  prob=transitions[hstates[i-1], ])
        seq[i] <- sample(emissions, 1, prob=emission_probs[hstates[i-1], ])
    }

    if (!is.null(dimnames(emission_probs))) {
        # make results labelled factor if emissoin probability matrix
        # has labels
        seq <- factor(seq, labels=colnames(emission_probs))
        hstates <- factor(hstates, labels=rownames(emission_probs))
    }
    list(sequence=seq, states=hstates)
}

### Coin Flip Simulation
## simulated biased coin
swaps_trans <- matrix(c(0.98, 0.02, 0.02, 0.98), nrow=2, ncol=2)
# H, T
initial_probs <- c(0.1, 0.9)
# (fair, biased) x (p(H), p(T))
coin_probs <- matrix(c(0.5, 0.8, 0.5, 0.2), nrow=2, ncol=2,
                     dimnames=list(c("F", "B"), c("H", "T")))

flips <- simulateSequence(100, swaps_trans, coin_probs, initial_probs)
