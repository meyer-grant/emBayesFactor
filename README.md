[![License](https://img.shields.io/badge/license-GPL(>=3)-C11B17.svg)](https://www.gnu.org/licenses/gpl-3.0.de.html)


# emBayesFactor
Computing Effect-Size and Moment Bayes Factors


## Description
This package provides functions for computing Bayes factors 
with effect-size and moment priors for one- and two-sample t tests, 
(multiple) linear regression models and analysis of variance (ANOVA). 


## Installation

### Using `devtools`
Make sure to install the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package
and then execute the command `devtools::install_github("meyer-grant/emBayesFactor@master")` in *R* 


## Examples
For examples please check out the help files for the `emBayesFactor` package in *R* via the commands
```
?ettest
?mttest
?ereg
?mreg
?eANOVA
?mANOVA
```


## Citation
If you want to cite this package write:

Meyer-Grant, C. G. & Klauer, K. C. (2024). *emBayesFactor: Computing effect-size and moment Bayes factors* [R package]. [https://github.com/meyer-grant/emBayesFactor](https://github.com/meyer-grant/emBayesFactor)


## References
Klauer, K. C., Meyer-Grant, C. G., & Kellen, D (in press). On Bayes Factors for Hypotheses Tests. *Psychonomic Bulletin & Review*. [PsyArXiv preprint: [https://doi.org/10.31234/osf.io/ykp29](https://doi.org/10.31234/osf.io/ykp29)]

Gronau, Q. F., Ly, A., & Wagenmakers, E.-J. (2020).Informed Bayesian *t*-tests. *The American Statistician*, *74*(2), 137–143.

Pramanik, S., & Johnson, V. E. (2024). Efficient alternatives for Bayesian hypothesis tests in psychology. *Psychological Methods*, *29*(2), 243–261.

Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. (2009). Bayesian *t*-tests for accepting and rejecting the null hypothesis. *Psychonomic Bulletin & Review*, *16*(2), 225–237.

