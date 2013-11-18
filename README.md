# hmmr -- HMMs in R

An implementation of the forward, backward, and foward-backward
algorithms in R, with scaling for numerical stability. While stable
and consistent with other library's results, you'll probably want to
use HMM in CRAN or Heng Li's
[khmm.c](https://github.com/attractivechaos/klib/blob/master/khmm.c)
for real data. These R implementations will eventually be wrapped in
calls to C.
