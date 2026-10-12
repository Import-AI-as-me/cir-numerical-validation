# CIR Pricing and Sensitivity Validation

**A price can match while its sensitivity is wrong.** In the analytic example below, the continuous-time CIR bond price and an annual-grid approximation both equal **0.8198175474**, yet the grid's initial-short-rate sensitivity is **3.06237% too high**.

This repository is a reproducible numerical validation study of CIR discounting. It separates integration bias, state-discretization bias, and Monte Carlo sampling error, with saved results and implementation checks. It reproduces established numerical mechanisms; it makes no claim of a new theorem or trading profitability.

[中文说明](README.zh-CN.md) · [English technical report](report.md) · [中文研究报告 PDF](report.zh-CN.pdf)

![An annual grid matches the price at a constructed initial rate while the initial-rate sensitivity remains about 3.06% too high.](price_greek_cancellation.png)

## The example

For the CIR short-rate model,

$$
dr_t=\kappa(\theta-r_t)\,dt+\sigma\sqrt{r_t}\,dW_t,
\qquad P=\mathbb E\!\left[e^{-\int_0^T r_s\,ds}\right],
\qquad G=-\frac{\partial P}{\partial r_0}.
$$

Use $\kappa=0.6$, $\theta=0.04$, $\sigma=0.1$, $T=5$ years, and an annual grid. At the deliberately constructed initial rate $r_0^*=0.0400853469$:

| Output | Continuous-time target | Exact grid transitions + trapezoidal integration |
| --- | ---: | ---: |
| Unit-notional price $P$ | 0.8198175474 | 0.8198175474 |
| Initial-rate sensitivity $G$ | 1.2852865306 | 1.3246467181 |
| Relative sensitivity error | — | +3.06237% |

The equality follows from cancellation of the affine pricing coefficients; the saved double-precision calculation reports zero absolute price difference. This is a constructed validation case, not an estimate of how often cancellation occurs in markets. Exact simulation of grid-point states still leaves error in the path integral. See [the saved counterexamples](cancellation_counterexamples.csv) and [the derivation](report.md).

## Run locally

Use **Python 3.12**, NumPy, and Matplotlib. All inputs are synthetic; no market-data subscription, API key, or network access is needed after dependency installation.

```bash
git clone https://github.com/Import-AI-as-me/cir-numerical-validation.git
cd cir-numerical-validation
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run the implementation tests and numerical checks (`--verify-only` prints checks as JSON without writing result files):

```bash
python -m unittest -v test_validation
python study.py --verify-only
```

The [recorded validation](validation.json) passed all 9 implementation tests, numerical checks, and the 36-configuration quick run in a fresh Python 3.12 environment.

Run a smaller experiment to inspect the workflow:

```bash
python study.py --quick --output-dir results/quick
```

Reproduce the full experiment configuration:

```bash
python study.py --reps 80 --paths 8000 --output-dir results/full
```

The quick run uses **4 replications of 512 paths** for each of the same 36 configurations. It is a workflow check and does **not** reproduce the full-study Monte Carlo tables. Runs refuse a nonempty output directory; choose a fresh directory for each rerun to preserve the committed reference results. [requirements-recorded.txt](requirements-recorded.txt) records the original environment; random-number streams, numerical rounding, and timings can differ across versions and machines.

## What is compared

All three schemes use the same trapezoidal approximation to the rate integral:

- **Exact:** noncentral chi-square CIR transitions at the chosen grid points.
- **Full truncation:** Euler updates with the positive part used in drift, diffusion, and discounting, while retaining the auxiliary state's negative values.
- **QE:** quadratic-exponential transitions matching the conditional mean and variance.

The main study uses three synthetic parameter sets, four grids (2, 4, 12, and 48 steps per year), and a five-year horizon. Each configuration has **80 independent replications of 8,000 paths**: 2,880 batches and 23.04 million paths overall. The base seed is `20261011`. The annual-grid cancellation example is a separate analytic calculation.

The affine formula for exact grid transitions provides the integration-bias reference. The Monte Carlo experiment measures **prices and price-interval coverage**; the sensitivity example compares analytic continuous-time and grid sensitivities. It does not validate Monte Carlo Greek estimators for Euler or QE.

## A result worth checking

In the synthetic `boundary` scenario, with a half-year step:

| Scheme | Estimated price error, unit-notional bp | MCSE of estimated error, bp | Coverage of the continuous-time price |
| --- | ---: | ---: | ---: |
| Exact | +1.397 | 1.114 | 76/80 |
| QE | +2.386 | 1.095 | 78/80 |
| Full truncation | −184.408 | 1.324 | 0/80 |

The analytic integration bias is only **+1.219 bp** here. The full-truncation intervals have a mean half-width of **23.546 bp**, so sampling error bars alone miss the dominant state-discretization error. More paths at a fixed biased grid do not remove that bias. Full data are in [summary.csv](summary.csv), with the broader comparison in [error_audit.png](error_audit.png).

Coverage is itself estimated. At 95% true coverage, 80 replications give an MCSE of about **2.44 percentage points**. Small differences between coverage percentages do not establish a reliable ranking. Observing 0/80 does not prove zero population coverage; its Wilson 95% upper bound is about 4.58%.

## Units and interpretation

- **Price-error bp** means $10^4(\widehat P-P)$ for unit notional. These are **not yield basis points**.
- **$G$** is sensitivity to the initial state $r_0$, holding $\kappa$, $\theta$, and $\sigma$ fixed. It is **not parallel-yield-curve DV01**. Legacy CSV columns named `*_r0_dv01` store $10^{-4}G$, the first-order price-decrease magnitude for a one-basis-point increase in $r_0$.
- **Bias MCSE** measures uncertainty in the estimated average error across independent batches. It is distinct from a single batch's price confidence interval.
- **Coverage** tests the usual mean ± 1.96 × estimated standard error interval against the continuous-time analytic price. It is not a finite-sample coverage guarantee.

The low-Feller cases are synthetic boundary stress tests, not calibrated market parameter estimates. Scheme comparisons use equal path counts, not equal runtime budgets. The study does not address parameter-estimation uncertainty, market calibration, or model misspecification.

## Repository contents

| File | Purpose |
| --- | --- |
| [study.py](study.py) | Simulation, affine references, checks, tables, and figures |
| [test_validation.py](test_validation.py) | Automated implementation checks |
| [report.md](report.md) / [report.zh-CN.pdf](report.zh-CN.pdf) | Derivations, experiment details, and limitations |
| [summary.csv](summary.csv) / [replicates.csv](replicates.csv) | Aggregate results / individual batch estimates |
| [analytic_grid.csv](analytic_grid.csv) | Exact grid pricing, sensitivities, and integration bias |
| [cancellation_counterexamples.csv](cancellation_counterexamples.csv) | Constructed price/sensitivity cancellation examples |
| [results.json](results.json) | Parameters, seed, environment, checks, and reference results |
| [requirements.txt](requirements.txt) / [requirements-recorded.txt](requirements-recorded.txt) | Installation requirements / recorded original environment |

## Relation to prior work

The grid-integral recursion is a transparent rederivation of existing work, especially **Peng and Kou (2008), §2.5.2, Lemma 2**. The confidence-interval distortion caused by discretization is also established; **Detemple, Garcia, and Rindisbacher (2003)** is a direct predecessor. This repository's purpose is to make these mechanisms inspectable in one reproducible pricing-and-sensitivity validation exercise.

- Cox, Ingersoll, and Ross (1985), [A Theory of the Term Structure of Interest Rates](https://doi.org/10.2307/1911242).
- Peng and Kou (2008), [Connecting the Top-down to the Bottom-up: Pricing CDO under a Conditional Survival (CS) Model](https://english.phbs.pku.edu.cn/uploadfile/2019/0920/20190920023124792.pdf).
- Detemple, Garcia, and Rindisbacher (2003), [Asymptotic Properties of Monte Carlo Estimators of Diffusion Processes](https://cirano.qc.ca/files/publications/2003s-11.pdf).

See the [technical report](report.md) for the numerical-scheme references and further research directions.

Development assistance: the initial draft and code were assisted by AI; reproducibility checks are included.
