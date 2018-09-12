data {
    int<lower=1> N;
    int<lower=0> n[N];
    vector<lower=0>[N] C; // offsets
}
parameters {
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
    n ~ poisson(adj_lambda);
}
