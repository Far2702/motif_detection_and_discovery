# **Statistical Motif Discovery in Regulatory DNA Sequences**

## **Overview**

This project implements a statistically principled framework for discovering and validating regulatory DNA motifs using classical bioinformatics models and hypothesis testing.
The approach constructs Position Frequency (PFM), Probability (PPM), and Weight Matrices (PWM) from curated regulatory sequences and scans input DNA sequences using log-likelihood ratio–based scoring.

To reduce false positives arising from dense, overlapping window scans, motif significance is evaluated under an empirical background model with conservative statistical thresholds.

An interactive Flask-based interface allows users to input sequences and explore statistically significant motif occurrences.


## Key Features

- Construction of PFM, PPM, and PWM from curated motif instances

- Log-likelihood ratio–based PWM scanning using sliding 22-mer windows

- Empirical background modeling for motif score distributions

- Hypothesis testing with α = 0.01 to control false positives

- Analysis of limitations of i.i.d. background assumptions

- Consideration of correlated hypothesis tests due to overlapping windows

- Interactive Flask-based web interface for motif scanning and visualization


## Methodology

- **Motif Modeling**

    - Construct PFM from aligned regulatory sequences

    - Convert PFM → PPM → PWM using background nucleotide frequencies

- **Motif Scanning**

    - Apply sliding-window scanning over query sequences

    - Compute log-likelihood ratio scores for each candidate site

- **Statistical Validation**

    - Estimate null distribution of scores using empirical background sequences

    - Perform hypothesis testing with conservative significance thresholds

    - Filter statistically significant motif occurrences

- **Analysis of Assumptions**

    - Examine independence assumptions in background models

    - Discuss effects of overlapping windows on multiple hypothesis testing

- **Planned Extensions**

    - Expectation–Maximization (EM)–based motif discovery from unlabeled sequences

    - Iterative refinement of PWM parameters

    - Support for variable-length motifs

    - Improved background models (Markov-based nucleotide dependencies)


## Technologies Used

- Python

- NumPy, Pandas, SciPy

- Flask

- Matplotlib / Plotly

 - Statistical modeling and hypothesis testing

 - Logomaker