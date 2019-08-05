// This is a adapted from of bc_mpra_model.stan (the main stan underlying malacoda)
// This Stan code defines a similar model, but it doesn't make any use of the priors, *_[m/p]_[a/b]
// This is so that rstan::optimizing() can be used to fit the MLEs
data {
  int<lower=0> n_rna_samples;
  int<lower=0> n_dna_samples;
  int<lower=1> n_ref; // number of reference barcodes
  int<lower=1> n_alt; // number of alternate barcodes
  int<lower=0> ref_dna_counts[n_ref, n_dna_samples];
  int<lower=0> alt_dna_counts[n_alt, n_dna_samples];
  int<lower=0> ref_rna_counts[n_ref, n_rna_samples];
  int<lower=0> alt_rna_counts[n_alt, n_rna_samples];

  real<lower=0> rna_depths[n_rna_samples]; 
  real<lower=0> dna_depths[n_dna_samples]; 
} 
parameters {
  real<lower=0> dna_m_ref;
  real<lower=0> dna_m_alt;
  real<lower=0> dna_p;

  vector<lower=0>[2] rna_m; // rna mean
  vector<lower=0>[2] rna_p; // rna phi aka size
}
model {

  // so here the individual counts come from a barcode-specific distribution, 
  // so the gamma prior on dna_m_a/b above needs to be fit on the 
  // depth adjusted counts themselves, not the mean
  for (s in 1:n_dna_samples) {
    for (bc in 1:n_ref) {
      ref_dna_counts[bc,s] ~ neg_binomial_2(dna_m_ref * dna_depths[s], dna_p);
    }
    
    for (bc in 1:n_alt) {
      alt_dna_counts[bc,s] ~ neg_binomial_2(dna_m_alt * dna_depths[s], dna_p);
    }
  }
  
  for (s in 1:n_rna_samples) {
    for (bc in 1:n_ref) {
      ref_rna_counts[bc, s] ~ neg_binomial_2(rna_m[1] * rna_depths[s] * dna_m_ref, rna_p[1]);
    }

    for (bc in 1:n_alt) {
      alt_rna_counts[bc, s] ~ neg_binomial_2(rna_m[2] * rna_depths[s] * dna_m_alt, rna_p[2]);
    }
  }

}
generated quantities {
  real ref_act;
  real alt_act;
  real transcription_shift;
  ref_act = log(rna_m[1]);
  alt_act = log(rna_m[2]);
  transcription_shift = alt_act - ref_act;
}
