data {
    int<lower=1> N;
    int<lower=0> n[N];
    real<lower=0> mu_init;
    vector<lower=0>[N] C; // offsets
}
parameters {
    real<lower=0,upper=1> theta;
    real<lower=0> alpha[2];
    real<lower=0> lambda[N];
}
transformed parameters {
    real<lower=0> adj_lambda[N];
    for (i in 1:N)
        adj_lambda[i] = lambda[i] * C[i];
}
model {
    alpha ~ gamma(0.01, 0.01);
    lambda ~ gamma(alpha[1], alpha[2]);
    for (i in 1:N) {
        if(n[i] < 1) {
            target += log_sum_exp(bernoulli_lpmf(1|theta),
                                  bernoulli_lpmf(0|theta) + poisson_lpmf(n[i]|adj_lambda[i]));
        } else {
            target += bernoulli_lpmf(0|theta) + poisson_lpmf(n[i]|adj_lambda[i]);
        }
    }
}
