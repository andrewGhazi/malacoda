data {
  int<lower=0> N;
  vector[N] mles;
  vector[N] weights;
}
parameters {
  real alpha;
  real beta;
}
model {
  for (i in 1:N) {
    target += weights[i] * gamma_lpdf(mles[i] | alpha, beta);
  }
}
