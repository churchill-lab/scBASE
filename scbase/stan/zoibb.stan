data {
    int<lower=1> N;
    int<lower=0> n[N];
    int<lower=0> x[N];
}
parameters {
    simplex[3] pi;
    real<lower=7> a_mono;
    vector<lower=2>[2] alpha;
    vector<lower=0,upper=1>[N] theta;
}
model {
    vector[3] log_pi = log(pi);
    a_mono ~ cauchy(7, 2);
    alpha ~ cauchy(2, 2);
    theta ~ beta(alpha[1], alpha[2]);
    for (i in 1:N) {
        vector[3] lps = log_pi;
        lps[1] = lps[1] + beta_binomial_lpmf(x[i]|n[i], a_mono, 1);
        lps[2] = lps[2] + beta_binomial_lpmf(x[i]|n[i], 1, a_mono);
        lps[3] = lps[3] + binomial_lpmf(x[i]|n[i], theta[i]);
        target += log_sum_exp(lps);
    }
}
generated quantities {
    matrix[N,3] pi_z;
    vector[N] theta_adj;
    real log_sum_exp_log_pi_z_raw;
    for (i in 1:N) {
        vector[3] log_pi_z_raw = log(pi);
        log_pi_z_raw[1] = log_pi_z_raw[1] + beta_binomial_lpmf(x[i]|n[i], a_mono, 1);
        log_pi_z_raw[2] = log_pi_z_raw[2] + beta_binomial_lpmf(x[i]|n[i], 1, a_mono);
        log_pi_z_raw[3] = log_pi_z_raw[3] + binomial_lpmf(x[i]|n[i], theta[i]);
        log_sum_exp_log_pi_z_raw = log_sum_exp(log_pi_z_raw);
        for (j in 1:3)
            pi_z[i,j] = exp(log_pi_z_raw[j] - log_sum_exp_log_pi_z_raw);
        theta_adj[i] = theta_adj[i] + pi_z[i, 1] * a_mono / (a_mono+1)
        theta_adj[i] = theta_adj[i] + pi_z[i, 2] * 1 / (a_mono+1)
        theta_adj[i] = theta_adj[i] + pi_z[i, 3] * theta[i]
    }
}
