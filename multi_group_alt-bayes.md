# Probabilistic multi-group comparison of alpha diversity
Rasmus Hindström
2025-07-24

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
# Model with partial pooling
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
    Intercept               1.46      0.14     1.19     1.74 1.00     9187     6229
    sigma_Intercept        -0.44      0.17    -0.76    -0.11 1.00     8043     5790
    AgeElderly             -0.15      0.27    -0.67     0.37 1.00     8341     6213
    AgeMiddle_age          -0.57      0.19    -0.95    -0.19 1.00     8836     6698
    sigma_AgeElderly        0.33      0.25    -0.15     0.83 1.00     9327     6510
    sigma_AgeMiddle_age    -0.27      0.26    -0.76     0.26 1.00     8483     6358

    Further Distributional Parameters:
       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    nu    24.40     14.71     5.70    60.36 1.00     9088     5644

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
library(ggplot2)
library(patchwork)

draws <- as_draws_df(fit)
population <- c(draws$b_Intercept, draws$b_Intercept + draws$b_AgeMiddle_age,  draws$b_Intercept + draws$b_AgeElderly)
pop_mean <- mean(population)

plot_data <- data.frame(
    pop_mean = pop_mean,
    adult = draws$b_Intercept,
    elderly = draws$b_Intercept + draws$b_AgeElderly,
    middle_age = draws$b_Intercept + draws$b_AgeMiddle_age
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

# 3. Classical approachs to multi-group testing

## 3.1. ANOVA

ANOVA would be the closest classical alternative. Assumptions in ANOVA
are normality and equal variance.

``` r
summary(aov(shannon ~ Age, data = df))
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
pairwise.t.test(df$shannon, df$Age, p.adjust = "fdr")
```


        Pairwise comparisons using t tests with pooled SD 

    data:  df$shannon and df$Age 

               Adult Elderly
    Elderly    0.505 -      
    Middle_age 0.047 0.132  

    P value adjustment method: fdr 

``` r
# HSD
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
kruskal.test(shannon ~ Age, df)
```


        Kruskal-Wallis rank sum test

    data:  shannon by Age
    Kruskal-Wallis chi-squared = 7.7239, df = 2, p-value = 0.02103

We get a p-value \< 0.05, so there are ‘significant’ differences within
the groups.

``` r
library(dunn.test)

dunn.test(df$shannon, df$Age, method = "bh")
```

      Kruskal-Wallis rank sum test

    data: x and group
    Kruskal-Wallis chi-squared = 7.7239, df = 2, p-value = 0.02

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

Significant p-values, after adjustment, are reported for the comparison
between groups `Adult` and `Middle age`.

# 4. Conclusions

Both the probabilistic and classical approachs to multigroup testing
give the same result. The `Middle age` group is signaled out.
Contrasting the methods, it is clear that the probabilistic approach
gives richer inference, with the added benefit of intuitive
interpretation.
