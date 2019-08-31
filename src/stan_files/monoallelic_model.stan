// This model is essentially a cut-down version of the main malacoda model, bc_mpra_model.stan.
// It models the RNA data for a single allele in the same way but doesn't go on to compare the alleles.
// There are some references to the "reference" allele, but this can be used to model either (mono) reference or alternate alleles.
// I was just too lazy to change the variable names.
data {
  int<lower=0> n_rna_samples;
  int<lower=0> n_dna_samples;
  int<lower=1> n_ref; // number of reference barcodes
  int<lower=0> ref_dna_counts[n_ref, n_dna_samples];
  int<lower=0> ref_rna_counts[n_ref, n_rna_samples];

  real<lower=0> rna_depths[n_rna_samples];
  real<lower=0> dna_depths[n_dna_samples];

  real<lower=0> dna_m_a; // ref and alt DNA use the same prior, hence length = 1
  real<lower=0> dna_m_b; // dna_m are basically the barcode effects here

  real<lower=0> dna_p_a; // dna_p priors are estimated differently
  real<lower=0> dna_p_b; // I should rename these variables

  real<lower=0> rna_m_a;
  real<lower=0> rna_m_b;
  real<lower=0> rna_p_a;
  real<lower=0> rna_p_b;
}
parameters {
  vector<lower=0>[n_ref] dna_m_ref;
  real<lower=0> dna_p;

  real<lower=0> rna_m; // rna mean
  real<lower=0> rna_p; // rna phi aka size
}
model {

  for (bc in 1:n_ref) {
    dna_m_ref[bc] ~ gamma(dna_m_a, dna_m_b);
  }

  dna_p ~ gamma(dna_p_a, dna_p_b);

  // so here the individual counts come from a barcode-specific distribution,
  // so the gamma prior on dna_m_a/b above needs to be fit on the
  // depth adjusted counts themselves, not the mean
  for (s in 1:n_dna_samples) {
    for (bc in 1:n_ref) {
      ref_dna_counts[bc,s] ~ neg_binomial_2(dna_m_ref[bc] * dna_depths[s], dna_p);
    }
  }

  rna_m ~ gamma(rna_m_a, rna_m_b); // priors on negative binomial parameters
  rna_p ~ gamma(rna_p_a, rna_p_b);

  for (s in 1:n_rna_samples) {
    for (bc in 1:n_ref) {
      ref_rna_counts[bc, s] ~ neg_binomial_2(rna_m * rna_depths[s] * dna_m_ref[bc], rna_p);
    }
  }

}
generated quantities {
  real activity;
  activity = log(rna_m);
}
