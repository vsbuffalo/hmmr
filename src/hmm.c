/*
A C HMM library by Vince Buffalo, <vsbuffaloAAAAA@gmail.com> sans-poly
A tail.

All the notation used follows Durbin et al, 1998 and Rabiner 1989.

Notation:
 - n states
 - m emitted symbols
 - state transition matrix A, where each element is a_{ij}
 - emissions probabilities, B={b_j(k)}
 - initial state probabilities

*/

#include "hmm.h"

hmm_t *init_hmm(int m, int n) {
  hmm_t hmm = malloc(sizeof(hmm_t));
  hmm->m = m;
  hmm->n = n;
  hmm->A = (double**) calloc(n*n, sizeof(double));
  hmm->B = (double**) calloc(n*m, sizeof(double));
  hmm->ip = calloc(n, sizeof(double));
  return hmm;
}

void destroy_hmm(hmm_t *hmm) {
  if (!hmm) return;
  free(hmm->A);
  hmm->B = (double**) calloc(n*m, sizeof(double));
  hmm->ip = calloc(n, sizeof(double));
  
}
