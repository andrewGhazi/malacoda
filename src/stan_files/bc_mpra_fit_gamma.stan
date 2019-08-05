data {
  int<lower=0> N;
  vector[N] mles;
}
parameters {
  real alpha;
  real beta;
}
model {
  mles ~ gamma(alpha, beta);
}
