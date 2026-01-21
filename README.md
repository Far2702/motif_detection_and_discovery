# **Statistical Motif Discovery and Detection in Regulatory DNA Sequences**

## Overview

This project implements a statistically principled framework for discovering and validating regulatory DNA motifs using classical bioinformatics models and hypothesis testing. The approach constructs Position Frequency (PFM), Position Probability (PPM), and Position Weight Matrices (PWM) from curated regulatory DNA sequences and scans query DNA sequences using log-likelihood ratio (LLR)–based scoring.

Motif significance is evaluated under an empirical background model to reduce false positives arising from dense, overlapping sliding-window scans. Conservative statistical thresholds are applied to ensure robust inference.

An interactive Flask-based interface allows users to input DNA sequences and explore statistically significant motif occurrences.

## Data Source

Motif instances were obtained from RegulonDB, a curated database of regulatory interactions.
Using this dataset:

- Point estimates of nucleotide proportions (A, T, G, C) were computed at each motif position

- These estimates were used to construct the PFM → PPM → PWM, forming the probabilistic motif model used throughout the analysis

## Key Features

- Construction of PFM, PPM, and PWM from curated motif instances

- Log-likelihood ratio–based PWM scanning using sliding 22-mer windows

- Empirical background modeling for motif score distributions

- Hypothesis testing with α = 0.01 to control false positives

- Analysis of limitations of i.i.d. background assumptions

- Consideration of correlated hypothesis tests due to overlapping windows

- Interactive Flask-based web interface for motif scanning and visualization

## Methodology

 ### Motif Modeling

- Construct PFM from aligned regulatory sequences

- Convert PFM → PPM → PWM using estimated background nucleotide frequencies

### Motif Scanning

- Apply sliding-window scanning across query DNA sequences

- Compute log-likelihood ratio (LLR) scores for each candidate subsequence

### Statistical Validation

- Evaluate motif significance using a hypothesis-testing framework

- Control false positives using conservative significance thresholds

### Analysis of Assumptions

- Examine independence assumptions in i.i.d. background models

- Discuss effects of overlapping windows on multiple hypothesis testing

## Statistical Hypothesis Testing Framework

### Test Statistic

The log-likelihood ratio (LLR) is used as the test statistic:

- Motif model: PWM-based probability

- Null model: Background nucleotide distribution

### Null Hypothesis (H₀)

<u>H₀:</u> *The observed LLR score of a 22-mer subsequence arises from the background nucleotide distribution alone and does not represent a true regulatory motif.*

### Alternative Hypothesis (H₁)

<u>H₁:</u> *The observed LLR score of a 22-mer subsequence is significantly higher than expected under the background model and is consistent with the PWM-defined regulatory motif.*

### Empirical Null Distribution Construction

Since the analytical distribution of the LLR random variable is unknown, statistical inference is performed using an empirical null distribution:

- Assume the null hypothesis (H₀) is true

- Generate many random DNA sequences from the background (null) model

- Compute LLR scores for each simulated sequence

- Aggregate these scores to obtain the empirical null distribution

This approach avoids parametric assumptions and allows direct estimation of significance thresholds.

### Decision Rule

- Choose a significance level α (default: α = 0.01)

- Reject H₀ if the observed LLR exceeds the (1 − α) quantile of the empirical null distribution

- Otherwise, fail to reject H₀

**Conservative α values are used due to correlated hypothesis tests induced by overlapping sliding windows.**

## Planned Extensions

- Expectation–Maximization (EM)–based motif discovery from unlabeled sequences

- Iterative refinement of PWM parameters

- Support for variable-length motifs

- Improved background models using Markov nucleotide dependencies


## Visualization and Interpretation

- **Sequence logos** derived from PPMs are used to visualize nucleotide conservation and degeneracy across motif positions

- **Interactive Plotly** visualizations are used to display empirical LLR distributions and significance thresholds, enabling intuitive interpretation of motif significance

## Technologies Used

- Python

- NumPy, Pandas, SciPy

- Flask

- Matplotlib / Plotly

- Logomaker

- Statistical modeling and hypothesis testing