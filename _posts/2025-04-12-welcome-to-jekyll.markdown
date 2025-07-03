---
title: Expectation Maximization Algorithm
layout: post
date: '2025-04-12 21:37:35 +0200'
categories: jekyll update
---

# Expectation-Maximization for ZINB Distributed scRNA-Seq Data

## Introduction

Single-cell RNA sequencing (scRNA-seq) generates sparse count matrices where 
many entries are zero. These zeros can result from both biological absence 
(True Negatives) and technical dropouts (False Negatives). To accurately model 
these features, the **Zero-Inflated Negative Binomial (ZINB)** model is often 
employed. 

This article explains how to use the Expectation-Maximization (EM) algorithm to 
estimate the parameters of the ZINB model in the context of scRNA-seq data 
analysis. 

### Objectives

- **Model Sparsity:** Handle excess zeros using a mixture model.
- **Parameter Estimation:** Use EM to estimate the zero-inflation probability 
alongside the Negative Binomial (NB) parameters (mean and dispersion).
- **Contextual Application:** Understand the application of mathematical methods 
within the framework of molecular biology and bioinformatics.

## The ZINB Model

In the ZINB model, each observed count y  is assumed to be generated from a 
two-component mixture:

1. **Zero-inflation component:** With probability \( \pi \), the observation is 
an "excess zero" (e.g., a dropout).
2. **Negative Binomial component:** With probability \( 1 - \pi \), the 
observation follows an NB distribution with parameters \( \mu \) (mean) and 
\( \theta \) (dispersion).

The probability mass function (pmf) is:

$$
P(y \mid \pi, \mu, \theta) = \begin{cases} \pi + (1 - \pi) \left( \dfrac{\theta}{\mu + \theta} \right)^{\theta}, & \text{if } y = 0, \\\\ (1 - \pi) \cdot \dfrac{\Gamma(y + \theta)}{y! \, \Gamma(\theta)} \left( \dfrac{\mu}{\mu + \theta} \right)^y \left( \dfrac{\theta}{\mu + \theta} \right)^{\theta}, & \text{if } y > 0. \end{cases}
$$

Here, \( \Gamma(\cdot) \) is the Gamma function.

## Formulating the EM Algorithm

Since we do not directly observe whether a zero is due to dropout or the NB process, we introduce a latent indicator variable \( z \):

- \( z = 1 \) if the count originates from the zero-inflation process.
- \( z = 0 \) if the count originates from the NB process.

### Complete-Data Log-Likelihood

For an observation \( y_i \) with latent indicator \( z_i \), the complete-data log-likelihood is given by:

$$
\log \mathcal{L}_c = \sum_{i} \left[ z_i \log \pi + (1 - z_i) \log (1 - \pi) + (1 - z_i) \log P_{\text{NB}}(y_i \mid \mu, \theta) \right],
$$

where \( P_{\text{NB}}(y_i \mid \mu, \theta) \) is the NB probability mass function.

## The EM Algorithm Steps

### Expectation (E) Step

For each observation \( y_i \), compute the expected latent indicator \( \gamma_i \), defined as follows:

- For \( y_i = 0 \):

  $$
  \gamma_i^{(t)} = \frac{\pi^{(t)}}{\pi^{(t)} + (1 - \pi^{(t)}) \left( \frac{\theta^{(t)}}{\mu^{(t)} + \theta^{(t)}} \right)^{\theta^{(t)}}}
  $$

- For \( y_i > 0 \):

  $$
  \gamma_i^{(t)} = 0,
  $$

since positive counts can only come from the NB component.

### Maximization (M) Step

Update the parameters by maximizing the expected complete-data log-likelihood:

1. **Update for \( \pi \):**

   $$
   \pi^{(t+1)} = \frac{\sum_{i : y_i = 0} \gamma_i^{(t)}}{N},
   $$

   where \( N \) is the total number of observations.

2. **Update for \( \mu \) and \( \theta \):**

   Numerical optimization (e.g., Newton-Raphson) is used to maximize:

   $$
   \sum_{i} (1 - \gamma_i^{(t)}) \log P_{\text{NB}}(y_i \mid \mu, \theta).
   $$

## Pseudocode Example

Below is a simplified pseudocode snippet (R-style) to illustrate the EM process:

```r
# Initialization
pi <- 0.5
mu <- mean(counts)
theta <- 1.0
tolerance <- 1e-6
converged <- FALSE

while (!converged) {
  # E-step: Calculate gamma for each observation
  gamma <- numeric(length(counts))
  for (i in 1:length(counts)) {
    if (counts[i] == 0) {
      nb_zero_prob <- (theta / (mu + theta))^theta
      gamma[i] <- pi / (pi + (1 - pi) * nb_zero_prob)
    } else {
      gamma[i] <- 0
    }
  }
  
  # M-step: Update parameter estimates
  pi_new <- sum(gamma) / length(counts)
  
  # Update mu and theta via numerical optimization (placeholder function)
  # [mu_new, theta_new] <- optimize_nb_parameters(counts, weights = (1 - gamma))
  
  # Check convergence
  if (abs(pi_new - pi) < tolerance) { # Apply similar checks for mu and theta
    converged <- TRUE
  }
  
  # Prepare for next iteration
  pi <- pi_new
  mu <- mu_new
  theta <- theta_new
}
```
