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

  real<lower=0> dna_m_a; // ref and alt DNA use the same prior, hence length = 1
  real<lower=0> dna_m_b; // dna_m are basically the barcode effects here

  real<lower=0> dna_p_a; // dna_p priors are estimated differently
  real<lower=0> dna_p_b; // I should rename these variables

  real<lower=0> rna_m_a[2]; // 1 = ref, 2 = alt
  real<lower=0> rna_m_b[2];
  real<lower=0> rna_p_a[2];
  real<lower=0> rna_p_b[2];
}
parameters {
  vector<lower=0>[n_ref] dna_m_ref;
  vector<lower=0>[n_alt] dna_m_alt;
  real<lower=0> dna_p;

  vector<lower=0>[2] rna_m; // rna mean
  vector<lower=0>[2] rna_p; // rna phi aka size
}
model {

  for (bc in 1:n_ref) {
    dna_m_ref[bc] ~ gamma(dna_m_a, dna_m_b);
  }

  for (bc in 1:n_alt) {
    dna_m_alt[bc] ~ gamma(dna_m_a, dna_m_b);
  }

  dna_p ~ gamma(dna_p_a, dna_p_b);

  // so here the individual counts come from a barcode-specific distribution,
  // so the gamma prior on dna_m_a/b above needs to be fit on the
  // depth adjusted counts themselves, not the mean
  for (s in 1:n_dna_samples) {
    for (bc in 1:n_ref) {
      ref_dna_counts[bc,s] ~ neg_binomial_2(dna_m_ref[bc] * dna_depths[s], dna_p);
    }

    for (bc in 1:n_alt) {
      alt_dna_counts[bc,s] ~ neg_binomial_2(dna_m_alt[bc] * dna_depths[s], dna_p);
    }
  }

  for (allele in 1:2) {
    rna_m[allele] ~ gamma(rna_m_a[allele], rna_m_b[allele]); // priors on negative binomial parameters
    rna_p[allele] ~ gamma(rna_p_a[allele], rna_p_b[allele]); // here, alleles have separate priors
  }

  for (s in 1:n_rna_samples) {
    for (bc in 1:n_ref) {
      ref_rna_counts[bc, s] ~ neg_binomial_2(rna_m[1] * rna_depths[s] * dna_m_ref[bc], rna_p[1]);
    }

    for (bc in 1:n_alt) {
      alt_rna_counts[bc, s] ~ neg_binomial_2(rna_m[2] * rna_depths[s] * dna_m_alt[bc], rna_p[2]);
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
