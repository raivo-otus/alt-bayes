# Probabilistic multi-group comparison of alpha diversity
Rasmus Hindström
2025-07-28

- [0. Summary](#0-summary)
- [1. Data preparation](#1-data-preparation)
- [2. Model fitting](#2-model-fitting)
- [3. Classical approachs to multi-group
  testing](#3-classical-approachs-to-multi-group-testing)
  - [3.1. ANOVA](#31-anova)
  - [3.2. Kruskal-Wallis](#32-kruskal-wallis)
- [4. Conclusions](#4-conclusions)

# 0. Summary

This report demonstrates the use and interpretation of a probabilisitc
alternative to multi-group comparisons of alpha diversity. The basis is
a multi-level model with a group-level effect, which is estimated using
Bayesian estimation with the `brms` package.

# 1. Data preparation

Preparing data from the `mia` package.

``` r
library(mia)
library(dplyr)
library(brms)
library(bayesplot)
library(dunn.test)
library(ggplot2)
library(patchwork)
```

``` r
data("peerj13075", package = "mia")
tse <- peerj13075
tse <- addAlpha(
    tse,
    assay.type = "counts",
    index = "shannon"
)
df <- as.data.frame(colData(tse))
```

# 2. Model fitting

The model is fitted using the `brm` function from the `brms` package.
Reponse variable is the Shannon index, and the grouping variable is the
3 classes of `age`; `Adult`, `Middle_age`, and `Elderly`.

The model is parametrized with the group `Adult` as the baseline to
which others are compared to.

Model definition is as follows:

$$
y_{ik} \sim \text{t}(\nu, \mu_{ik}, \sigma_k)
$$

$$
\mu_{ik} = \beta_0 + \beta_k
$$

$$
\sigma_k = \gamma_0 + \gamma_k
$$

*Default priors used by `brm()`*

$$
\nu \sim \gamma(2, 0.1)
$$ $$
\beta_0 \sim \text{t}(3, 1.3, 2.5)
$$ $$
\gamma_0 \sim \text{t}(3, 0, 2.5)
$$

$$
\beta_k, \gamma_k \sim \mathrm{Uniform}
$$

``` r
start <- proc.time()
fit <- brm(
    formula = bf(
        shannon ~ Age,
        sigma ~ Age
    ),
    data = df,
    family = student(),
    iter = 4000,
    chains = 4,
    cores = 4
)
end <- proc.time()
runTime_brm <- end - start
```

     Family: student 
      Links: mu = identity; sigma = log; nu = identity 
    Formula: shannon ~ Age 
             sigma ~ Age
       Data: df (Number of observations: 58) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    Intercept               1.46      0.14     1.19     1.73 1.00     8577     5755
    sigma_Intercept        -0.44      0.17    -0.76    -0.10 1.00     7054     5534
    AgeElderly             -0.14      0.26    -0.67     0.38 1.00     7780     6190
    AgeMiddle_age          -0.57      0.19    -0.95    -0.19 1.00     7929     6383
    sigma_AgeElderly        0.33      0.25    -0.15     0.84 1.00     7629     6257
    sigma_AgeMiddle_age    -0.26      0.26    -0.76     0.24 1.00     8198     5482

    Further Distributional Parameters:
       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    nu    24.24     14.34     5.87    59.77 1.00     8680     5134

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

Interpreting the coefficients and 95% CI we can already make the
observation that the groups `Adult` and `Middle age` differ. The
`Middle age` group appears to have a lower Shannon diversity.

Further plotting is required to make conclusions on other pair wise
comparisons.

<details class="code-fold">
<summary>Posterior plotting</summary>

``` r
draws <- as_draws_df(fit)
population <- c(draws$b_Intercept, draws$b_Intercept + draws$b_AgeMiddle_age,  draws$b_Intercept + draws$b_AgeElderly)
pop_mean <- mean(population)

plot_data <- data.frame(
    pop_mean = pop_mean,
    adult = draws$b_Intercept,
    adult_sd = draws$b_sigma_Intercept,
    elderly = draws$b_Intercept + draws$b_AgeElderly,
    elderly_sd = draws$b_sigma_Intercept + draws$b_sigma_AgeElderly,
    middle_age = draws$b_Intercept + draws$b_AgeMiddle_age,
    middle_age_sd = draws$b_sigma_Intercept + draws$b_sigma_AgeMiddle_age
)

p1 <- ggplot(data = plot_data) +
    geom_density(aes(x = adult), fill = "blue", alpha = 0.5, color = "blue") +
    geom_density(aes(x = middle_age), fill = "orange", alpha = 0.6, color = "orange") +
    geom_density(aes(x = elderly), fill = "purple", alpha = 0.7, color = "purple") +
    geom_vline(xintercept = plot_data$pop_mean, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
        title = "Posterior Distributions of Group means",
        x = "Shannon index",
        y = "Density"
    ) +
    annotate(
        "text", x = 2, y = 2.5,
        label = "Blue: Adult\nOrange: Middle age\nPurple: Elderly"
    ) +
    theme_minimal()

p2 <- ggplot(data = plot_data) +
    geom_boxplot(aes(y = adult, x = 1), fill = "blue", alpha = 0.7, color = "black") +
    geom_boxplot(aes(y = middle_age, x = 2), fill = "orange", alpha = 0.7, color = "black") +
    geom_boxplot(aes(y = elderly, x = 3), fill = "purple", alpha = 0.7, color = "black") +
    geom_hline(yintercept = plot_data$pop_mean, linetype = "dashed", color = "red", linewidth = 1) +
    scale_x_continuous(
        breaks = c(1, 2, 3),
        labels = c("Adult", "Middle Age", "Elderly")
    ) +
    labs(
        title = "Boxplots of Group Means",
        x = "Group",
        y = "Shannon Index"
    ) +
    theme_minimal()

p1 + p2
```

</details>

![](multi_group_alt-bayes_files/figure-commonmark/plotting-post-1.png)

From the plots we can infer groups are not similar. Particularly the
Middle aged (Orange) group appears to have a lower Shannon index. The
boxplots paint a clear picture of the higher overlap between the Adult
and Edlerly group, while the Middle aged group differs. In both plots
the red dashed line indicates the total population posterior mean.

Using the posterior distributions, we can make statements about the
differences between the groups. In this context the probability of
observing a higher Shannon index is appropriate, and akin to a
frequentist p-value.

<details class="code-fold">
<summary>Probabilities and Standardized Effect Size</summary>

``` r
calc_eff_posterior <- function(mu_1, mu_2, sd_1, sd_2) {
    # Calculates effect size from posterior draws
    diff <- mu_1 - mu_2
    pooled_sd <- sqrt((sd_1^2 + sd_2^2) / 2)
    d <- diff / pooled_sd
    mean_d <- mean(d)
    ci_d <- quantile(d, probs = c(0.05, 0.95))
    res <- list(
        "d" = mean_d,
        "ci" = ci_d
    )
    return(res)
}

probabilities <- data.frame(
    Comparison = c(
        "Adult vs Elderly",
        "Adult vs Middle age",
        "Elderly vs Middle age"
    ),
    Prob_greater = c(
        prob_adult_elderly <- mean(plot_data$adult > plot_data$elderly),
        prob_adult_middleage <- mean(plot_data$adult > plot_data$middle_age),
        prob_elderly_middleage <- mean(plot_data$elderly > plot_data$middle_age)
    ),
    High_P = c(
        ifelse(prob_adult_elderly > 0.95, "*", ""),
        ifelse(prob_adult_middleage > 0.95, "*", ""),
        ifelse(prob_elderly_middleage > 0.95, "*", "")
    ),
    cohens_d = c(
        calc_eff_posterior(plot_data$adult, plot_data$elderly, plot_data$adult_sd, plot_data$elderly_sd)$d,
        calc_eff_posterior(plot_data$adult, plot_data$middle_age, plot_data$adult_sd, plot_data$middle_age_sd)$d,
        calc_eff_posterior(plot_data$elderly, plot_data$middle_age, plot_data$elderly_sd, plot_data$middle_age_sd)$d
    ),
    ci_lower = c(
        calc_eff_posterior(plot_data$adult, plot_data$elderly, plot_data$adult_sd, plot_data$elderly_sd)$ci[1],
        calc_eff_posterior(plot_data$adult, plot_data$middle_age, plot_data$adult_sd, plot_data$middle_age_sd)$ci[1],
        calc_eff_posterior(plot_data$elderly, plot_data$middle_age, plot_data$elderly_sd, plot_data$middle_age_sd)$ci[1]
    ),
    ci_upper = c(
        calc_eff_posterior(plot_data$adult, plot_data$elderly, plot_data$adult_sd, plot_data$elderly_sd)$ci[2],
        calc_eff_posterior(plot_data$adult, plot_data$middle_age, plot_data$adult_sd, plot_data$middle_age_sd)$ci[2],
        calc_eff_posterior(plot_data$elderly, plot_data$middle_age, plot_data$elderly_sd, plot_data$middle_age_sd)$ci[2]
    )
)

knitr::kable(probabilities, caption = "Probabilities of Higher Shannon Index", format = "pipe")
```

</details>

| Comparison            | Prob_greater | High_P |  cohens_d |   ci_lower | ci_upper |
|:----------------------|-------------:|:-------|----------:|-----------:|---------:|
| Adult vs Elderly      |     0.710750 |        | 0.4852815 | -0.9502119 | 2.133642 |
| Adult vs Middle age   |     0.997500 | \*     | 0.9909632 |  0.4073328 | 1.788763 |
| Elderly vs Middle age |     0.944625 |        | 0.8840643 | -0.0258875 | 2.054133 |

Probabilities of Higher Shannon Index

Notice, that these probabilities are not p-values, and their
interpretation is more intuitive. The probabilities can also be computed
over the highest density interval’s (HDI) or 95% CI’s. Here we have used
the full posterior distributions.

Effect size’s are standardized and 95% CI’s are reported.

# 3. Classical approachs to multi-group testing

## 3.1. ANOVA

ANOVA would be the closest classical alternative. Assumptions in ANOVA
are normality and equal variance.

``` r
start <- proc.time() 
res_annova <- aov(shannon ~ Age, data = df)
end <- proc.time()
runTime_annova <- end - start

summary(res_annova)
```

                Df Sum Sq Mean Sq F value Pr(>F)  
    Age          2  3.188  1.5942   3.194 0.0487 *
    Residuals   55 27.448  0.4991                 
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Our result inidcates that Age is a significant predictor on explaining
differences in the Shannon index. After observing a significant result,
it is typical to further inspect the data with a pairwise t-test.
Another option is Tukey’s Honest Significant Difference (HSD). Both of
which require adjusting p-values due to being post-hoc tests.

``` r
# Pairwise t-test
start <- proc.time()
pairwise.t.test(df$shannon, df$Age, p.adjust = "fdr")
```


        Pairwise comparisons using t tests with pooled SD 

    data:  df$shannon and df$Age 

               Adult Elderly
    Elderly    0.505 -      
    Middle_age 0.047 0.132  

    P value adjustment method: fdr 

``` r
end <- proc.time()
runTime_pwt <- end - start

# HSD
start <- proc.time()
TukeyHSD(aov(shannon ~ Age, data = df))
```

      Tukey multiple comparisons of means
        95% family-wise confidence level

    Fit: aov(formula = shannon ~ Age, data = df)

    $Age
                             diff        lwr         upr     p adj
    Elderly-Adult      -0.1477303 -0.6783059  0.38284533 0.7814198
    Middle_age-Adult   -0.5690925 -1.1182905 -0.01989461 0.0406666
    Middle_age-Elderly -0.4213623 -1.0060281  0.16330359 0.2011258

``` r
end <- proc.time()
runTime_hsd <- end - start
```

Both post-hoc tests point to the significant difference between the
groups `Adult` and `Middle age`. In the pairwise t-test the difference
between groups `Middle age` and `Elderly` approaches the significance
boundry (p \< 0.05).

## 3.2. Kruskal-Wallis

Another option to ANOVA is the Kruskal-Wallis test. Kruskal-Wallis
relaxes the assumption of normality. However, It also needs to be paired
with a post-hoc test to infer paired differences after global
significance has been tested. Kruskal-Wallis is typically paired with
Dunn’s post-hoc test.

``` r
start <- proc.time()
kruskal.test(shannon ~ Age, df)
```


        Kruskal-Wallis rank sum test

    data:  shannon by Age
    Kruskal-Wallis chi-squared = 7.7239, df = 2, p-value = 0.02103

``` r
end <- proc.time()
runTime_kw <- end - start
```

We get a p-value \< 0.05, so there are ‘significant’ differences within
the groups.

``` r
start <- proc.time()
dunn.test(df$shannon, df$Age, method = "bh", kw = FALSE)
```


                               Comparison of x by group                            
                                 (Benjamini-Hochberg)                              
    Col Mean-|
    Row Mean |      Adult    Elderly
    ---------+----------------------
     Elderly |   0.672628
             |     0.2506
             |
    Middle_a |   2.733071   1.956873
             |    0.0094*     0.0378

    alpha = 0.05
    Reject Ho if p <= alpha/2

``` r
end <- proc.time()
runTime_dunn <- end - start
```

Significant p-values, after adjustment, are reported for the comparison
between groups `Adult` and `Middle age`.

# 4. Conclusions

<details class="code-fold">
<summary>Comparison of methods</summary>

``` r
runTimes <- data.frame(
    method = c(
        "Bayesian estimation",
        "ANNOVA + t.test",
        "ANNOVA + HSD",
        "Kruskal-Wallis + Dunn's"
        ),
    time_seconds = c(
        runTime_brm["elapsed"],
        runTime_annova["elapsed"] + runTime_pwt["elapsed"],
        runTime_annova["elapsed"] + runTime_hsd["elapsed"],
        runTime_kw["elapsed"] + runTime_dunn["elapsed"]
        )
)

knitr::kable(runTimes, caption = "Run times for different methods", format = "pipe")
```

</details>

| method                  | time_seconds |
|:------------------------|-------------:|
| Bayesian estimation     |       87.986 |
| ANNOVA + t.test         |        0.004 |
| ANNOVA + HSD            |        0.009 |
| Kruskal-Wallis + Dunn’s |        0.009 |

Run times for different methods

Both the probabilistic and classical approachs to multigroup testing
give the same result. The `Middle age` group is signaled out.
Contrasting the methods, it is clear that the probabilistic approach
gives richer inference, with the added benefit of intuitive
interpretation.

The main drawback of the probabilistic approach is the computational
cost. Fitting a Bayesian model, with brms, requires compiling the model
and running the MCMC sampler.
