// This is a tissue model of MPRA data.
// It's adapted from bc_mpra_model but allows for 1 OR 2 alleles and an arbitrary number of tissues.
// It DOESN'T account for an interaction between tissue and allele
data {
  int<lower=0> n_rna_samples; // the TOTAL numer of RNA samples
  int<lower=0> n_dna_samples;
  int<lower=1> n_bc;
  int<lower=1> n_allele; // Putting no upper limit -- I guess this can handle tri/quad allelic variants
  int<lower=1> n_tissue;
  int<lower=1> n_samples_per_tissue[n_tissue];

  int<lower=1, upper=n_bc> dna_bc_id[n_bc];
  int<lower=1, upper=n_bc> rna_bc_id[n_bc];
  int<lower=1, upper=n_allele> allele_id[n_bc];
  int<lower=1, upper=n_rna_samples> tissue_id[n_rna_samples];
  int<lower=1, upper=n_tissue> sample_tissue_id[n_rna_samples];

  int<lower=0> dna_counts[n_bc,n_dna_samples];
  int<lower=0> rna_counts[n_bc,n_rna_samples];

  real<lower=0> rna_depths[n_rna_samples];
  real<lower=0> dna_depths[n_dna_samples];

  real<lower=0> dna_m_a; // ref and alt DNA use the same prior, hence length = 1
  real<lower=0> dna_m_b; // dna_m are basically the barcode effects here

  real<lower=0> dna_p_a; // dna_p priors are estimated differently
  real<lower=0> dna_p_b; // I should rename these variables

  real<lower=0> rna_m_a[n_tissue,n_allele];
  real<lower=0> rna_m_b[n_tissue,n_allele];
  real<lower=0> rna_p_a[n_tissue,n_allele];
  real<lower=0> rna_p_b[n_tissue,n_allele];
}
parameters {
  vector<lower=0>[n_bc] dna_m;
  real<lower=0> dna_p;

  vector[n_tissue-1] tissue_effect; // the first tissue and allele are the "reference" class
  vector[n_allele-1] allele_effect;
  real<lower=0> rna_m; // rna mean
  real<lower=0> rna_p; // rna phi aka size
  // ^ Let's assume the RNA dispersions are all the same for now
  // And that all variation in RNA mean is due to depth, tissue, DNA, and allele
}
model {

  dna_m ~ gamma(dna_m_a, dna_m_b);
  dna_p ~ gamma(dna_p_a, dna_p_b);

  for (t in 1:(n_tissue-1)){
    tissue_effect[t] ~ student_t(5, 0, 1);
  }

  for (a in 1:(n_allele - 1)){
    allele_effect[a] ~ student_t(5, 0, 1);
  }

  for (t in 1:n_tissue){
    for (a in 1:n_allele){
      rna_m ~ gamma(rna_m_a[t,a], rna_m_b[t,a]);
      rna_p ~ gamma(rna_p_a[t,a], rna_p_b[t,a]);
    }
  }


  for (s in 1:n_dna_samples){
    for (bc in 1:n_bc){
      dna_counts[bc, s] ~ neg_binomial_2(dna_m[bc] * dna_depths[s], dna_p);
    }
  }

  for (s in 1:n_rna_samples){
    for (bc in 1:n_bc){
      if (allele_id[bc] == 1 && sample_tissue_id[s] == 1){
        rna_counts[bc, s] ~ neg_binomial_2(dna_m[bc] * rna_depths[s] * rna_m, rna_p);
      } else if (allele_id[bc] == 1 && sample_tissue_id[s] != 1) {
        rna_counts[bc, s] ~ neg_binomial_2(dna_m[bc] * rna_depths[s] * exp(tissue_effect[sample_tissue_id[s]-1]) * rna_m, rna_p);
      } else if (allele_id[bc] != 1 && sample_tissue_id[s] == 1) {
        rna_counts[bc, s] ~ neg_binomial_2(dna_m[bc] * rna_depths[s] * exp(allele_effect[allele_id[bc]-1]) * rna_m, rna_p);
      } else {
        rna_counts[bc, s] ~ neg_binomial_2(dna_m[bc] * rna_depths[s] * exp(tissue_effect[sample_tissue_id[s]-1]) * exp(allele_effect[allele_id[bc]-1]) * rna_m, rna_p);
      }
    }
  }
}

