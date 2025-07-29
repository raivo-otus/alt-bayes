# Probabilistic multi-group comparison of alpha diversity
Rasmus Hindström
2025-07-29

- [0. Summary](#0-summary)
- [1. Data preparation](#1-data-preparation)
- [2. Model fitting](#2-model-fitting)
- [3. Posterior plotting](#3-posterior-plotting)
- [4. Quantify probabilities and effect
  sizes](#4-quantify-probabilities-and-effect-sizes)
- [5. Classical approachs to multi-group
  testing](#5-classical-approachs-to-multi-group-testing)
  - [5.1. ANOVA](#51-anova)
  - [5.2. Kruskal-Wallis](#52-kruskal-wallis)
- [6. Conclusions](#6-conclusions)

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
    Intercept               1.46      0.14     1.18     1.73 1.00     9424     5770
    sigma_Intercept        -0.45      0.17    -0.76    -0.12 1.00     8332     5887
    AgeElderly             -0.14      0.27    -0.68     0.38 1.00     8612     5900
    AgeMiddle_age          -0.57      0.19    -0.94    -0.18 1.00     9629     5904
    sigma_AgeElderly        0.34      0.25    -0.13     0.83 1.00     8496     6407
    sigma_AgeMiddle_age    -0.26      0.25    -0.76     0.26 1.00     9142     6310

    Further Distributional Parameters:
       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    nu    24.33     14.47     5.88    60.72 1.00    10603     5646

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

Interpreting the coefficients and 95% CI we can already make the
observation that the groups `Adult` and `Middle age` differ. The
`Middle age` group appears to have a lower Shannon diversity.

Further plotting is required to make conclusions on other pair wise
comparisons.

# 3. Posterior plotting

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
observing a lower Shannon index is appropriate, akin to a classical
p-value.

# 4. Quantify probabilities and effect sizes

<details class="code-fold">
<summary>Probabilities and Standardized Effect Size</summary>

``` r
calc_cohenD_post <- function(mu_1, mu_2, sd_1, sd_2) {
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

calc_logfc_post <- function(mu_1, mu_2) {
    # Calculates log fold change from posterior draws
    logfc <- log2(mu_1 / mu_2)
    mean_logfc <- mean(logfc)
    ci_logfc <- quantile(logfc, probs = c(0.05, 0.95))
    res <- list(
        "logfc" = mean_logfc,
        "ci" = ci_logfc
    )
    return(res)
}

# Calculate effect sizes
logfc_adult_elderly <- calc_logfc_post(plot_data$adult, plot_data$elderly)
logfc_adult_middleage <- calc_logfc_post(plot_data$adult, plot_data$middle_age)
logfc_elderly_middleage <- calc_logfc_post(plot_data$elderly, plot_data$middle_age)

cohenD_adult_elderly <- calc_cohenD_post(
    plot_data$adult, 
    plot_data$elderly, 
    plot_data$adult_sd, 
    plot_data$elderly_sd
)
cohenD_adult_middleage <- calc_cohenD_post(
    plot_data$adult, 
    plot_data$middle_age, 
    plot_data$adult_sd, 
    plot_data$middle_age_sd
)
cohenD_elderly_middleage <- calc_cohenD_post(
    plot_data$elderly, 
    plot_data$middle_age, 
    plot_data$elderly_sd, 
    plot_data$middle_age_sd
)



probabilities <- data.frame(
    Comparison = c(
        "Adult vs Elderly",
        "Adult vs Middle age",
        "Elderly vs Middle age"
    ),
    Prob_lesser = c(
        prob_adult_elderly <- mean(plot_data$adult < plot_data$elderly),
        prob_adult_middleage <- mean(plot_data$adult < plot_data$middle_age),
        prob_elderly_middleage <- mean(plot_data$elderly < plot_data$middle_age)
    ),
    LogFC = c(
        logfc_adult_elderly$logfc,
        logfc_adult_middleage$logfc,
        logfc_elderly_middleage$logfc
    ),
    LogFC_ci_lower = c(
        logfc_adult_elderly$ci[1],
        logfc_adult_middleage$ci[1],
        logfc_elderly_middleage$ci[1]
    ),
    LogFC_ci_upper = c(
        logfc_adult_elderly$ci[2],
        logfc_adult_middleage$ci[2],
        logfc_elderly_middleage$ci[2]
    ),
    cohens_d = c(
        cohenD_adult_elderly$d,
        cohenD_adult_middleage$d,
        cohenD_elderly_middleage$d
    ),
    d_ci_lower = c(
        cohenD_adult_elderly$ci[1],
        cohenD_adult_middleage$ci[1],
        cohenD_elderly_middleage$ci[1]
    ),
    d_ci_upper = c(
        cohenD_adult_elderly$ci[2],
        cohenD_adult_middleage$ci[2],
        cohenD_elderly_middleage$ci[2]
    )
)

knitr::kable(probabilities, caption = "", format = "pipe")
```

</details>

| Comparison | Prob_lesser | LogFC | LogFC_ci_lower | LogFC_ci_upper | cohens_d | d_ci_lower | d_ci_upper |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Adult vs Elderly | 0.287875 | 0.1668741 | -0.2885604 | 0.6715804 | 0.4911642 | -0.9590143 | 2.136760 |
| Adult vs Middle age | 0.002625 | 0.7183928 | 0.3130507 | 1.1503023 | 0.9875781 | 0.4123324 | 1.795043 |
| Elderly vs Middle age | 0.051500 | 0.5515186 | -0.0119695 | 1.1022164 | 0.8849451 | -0.0145154 | 2.068362 |

Notice, that these are not classical p-values, but posterior
probabilities. The probabilities of observing a lower shannon index in
the first group of each comparison are reported. They have been
calculated from the full posterior distribution.

Effect size’s are reported as Log Fold Change (LogFC) and standardized
effect size (Cohen’s d) with 95% CI’s.

Combining the probabilities and effect sizes, we can conclude that the
`Middle age` group has a lower Shannon index compared to the `Adult`
group. With high probability (0.99) of observing a higher Shannon index
in the `Adult` group and a higher then 95% chance of observing an
positive effect size.

The `Elderly` group has a similar pattern, with a high probabilty of
observing a higher Shannon index than the `Middle age` group, but effect
size 95% CI’s overlap zero. This indicates that the difference is not as
pronounced as with the `Adult` group.

Comparing the `Adult` and `Elderly` groups, the probability of observing
a higher Shannon index in the `Adult` group is greater then 50%, but the
effect size is small, with a 95% CI that overlaps zero. This indicates
that the difference is not likely to be meaningful.

# 5. Classical approachs to multi-group testing

## 5.1. ANOVA

ANOVA would be the closest classical alternative. Assumptions in ANOVA
are normality and equal variance.

``` r
start <- proc.time() 
res_anova <- aov(shannon ~ Age, data = df)
end <- proc.time()
runTime_anova <- end - start

summary(res_anova)
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

## 5.2. Kruskal-Wallis

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

# 6. Conclusions

<details class="code-fold">
<summary>Comparison of run times</summary>

``` r
runTimes <- data.frame(
    method = c(
        "Bayesian estimation",
        "ANOVA + t.test",
        "ANOVA + HSD",
        "Kruskal-Wallis + Dunn's"
        ),
    time_seconds = c(
        runTime_brm["elapsed"],
        runTime_anova["elapsed"] + runTime_pwt["elapsed"],
        runTime_anova["elapsed"] + runTime_hsd["elapsed"],
        runTime_kw["elapsed"] + runTime_dunn["elapsed"]
        )
)

knitr::kable(runTimes, caption = "", format = "pipe")
```

</details>

| method                  | time_seconds |
|:------------------------|-------------:|
| Bayesian estimation     |       83.126 |
| ANOVA + t.test          |        0.005 |
| ANOVA + HSD             |        0.007 |
| Kruskal-Wallis + Dunn’s |        0.006 |

Although the methods are in aggreement that the significant outlier is
the `Middle age` group. The classical methods are only able to provide a
binary (yes/no) inference, in contrast to the much richer and intuitive
information provided by the probabilistic method.

The main drawback of the probabilistic approach is the computational
cost. Fitting a Bayesian model, with brms, requires compiling the model
code and running the MCMC sampler. However this margin may shrink when
datasets grow larger, as probabilistic model fitting scales well larger
data.
