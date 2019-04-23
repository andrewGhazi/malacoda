data {
  int<lower=0> n_grna; // number of gRNAs for the gene in the assay
  int<lower=0> n_in; // number input sequencing samples
  int<lower=0> n_out; // number output sequencing samples

  int<lower=0> in_counts[n_grna, n_in];
  int<lower=0> out_counts[n_grna, n_out];

  real<lower=0> in_depths[n_in];
  real<lower=0> out_depths[n_out];

  // The counts are modeled as coming from negative binomial distributions.
  // NB's have two parameters: mean and size.
  // Both the mean and size get their own gamma distribution priors.
  // The gamma distribution has two parameters: alpha and beta.
  // So eight inputs are required to specify the priors for this model.

  real<lower=0> in_mean_a; // input mean alpha
  real<lower=0> in_mean_b; // input mean beta
  real<lower=0> in_size_a; // input size alpha
  real<lower=0> in_size_b; // input size alpha

  real<lower=0> out_mean_a; // likewise
  real<lower=0> out_mean_b;
  real<lower=0> out_size_a;
  real<lower=0> out_size_b;
}

parameters {
  real<lower=0> in_mean[n_grna]; // the "concentration" of each input gRNA is modeled
  real<lower=0> in_size;
  real<lower=0> out_mean; // there's only one out_mean since they should be the same after accounting for depth & input concentration
  real<lower=0> out_size;
}

model {

  // gamma priors on NB parameters
  for (g in 1:n_grna) {
    in_mean[g] ~ gamma(in_mean_a, in_mean_b);
  }
  in_size ~ gamma(in_size_a, in_size_b);
  out_mean ~ gamma(out_mean_a, out_mean_b);
  out_size ~ gamma(out_size_a, out_size_b);

  // NB distributions on counts
  for (s in 1:n_in){
    for (g in 1:n_grna){
      in_counts[g,s] ~ neg_binomial_2(in_mean[g] * in_depths[s], in_size);
    }
  }

  for (s in 1:n_out){
    for (g in 1:n_grna){
      out_counts[g,s] ~ neg_binomial_2(in_mean[g] * out_depths[s] * out_mean, out_size);
    }
  }
}
