data {
  int<lower=1> n_gRNA;

  int<lower=0> t0_counts[n_gRNA];
  int<lower=0> case_counts[n_gRNA];
  int<lower=0> ctrl_counts[n_gRNA];

  real<lower=0> t0_depths;
  real<lower=0> case_depths;
  real<lower=0> ctrl_depths;

  real<lower=0> t0_m_a;
  real<lower=0> t0_m_b;

  real<lower=0> t0_p_a; // t0_p priors are estimated differently
  real<lower=0> t0_p_b; // I should rename these variables

  real<lower=0> output_m_a[2]; // 1 = ctrl, 2 = case
  real<lower=0> output_m_b[2];
  real<lower=0> output_p_a[2];
  real<lower=0> output_p_b[2];
}
parameters {
  vector<lower=0>[n_gRNA] t0_m;
  real<lower=0> t0_p;

  vector<lower=0>[2] output_m; // 1 = ctrl, 2 = case
  vector<lower=0>[2] output_p;

}
model {

  for (g in 1:n_gRNA) {
    t0_m[g] ~ gamma(t0_m_a, t0_m_b);
  }

  t0_p ~ gamma(t0_p_a, t0_p_b);

  // so here the individual counts come from a barcode-specific distribution,
  // so the gamma prior on dna_m_a/b above needs to be fit on the
  // depth adjusted counts themselves, not the mean
  for (g in 1:n_gRNA) {
    t0_counts[g] ~ neg_binomial_2(t0_m[g] * t0_depths, t0_p);
  }


  for (group in 1:2) {
    output_m[group] ~ gamma(output_m_a[group], output_m_b[group]); // priors on negative binomial parameters
    output_p[group] ~ gamma(output_p_a[group], output_p_b[group]); // here, groups have separate priors
  }

for (g in 1:n_gRNA) {
      ctrl_counts[g] ~ neg_binomial_2(output_m[1] * ctrl_depths * t0_m[g], output_p[1]);
    }

    for (g in 1:n_gRNA) {
      case_counts[g] ~ neg_binomial_2(output_m[2] * case_depths * t0_m[g], output_p[2]);
    }


}
generated quantities {
  real ctrl_act;
  real case_act;
  real condition_shift;
  ctrl_act = log(output_m[1]);
  case_act = log(output_m[2]);
  condition_shift = case_act - ctrl_act;
}
