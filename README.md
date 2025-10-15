# CFM-PLV: Bayesian circular functional models for computing noise-robust phase locking values

This repository provides R code implementing noise-robust phase locking values based on Bayesian circular functional models, as proposed by the following paper.

Sugasawa, S., Matsuda, T. and Nakagawa, T. (2025).  Noise-robust phase connectivity estimation via Bayesian circular functional models. [arXiv:2509.06418](https://arxiv.org/abs/2509.06418)

This repository includes the following files: 

- `WMF.R`: Scriot implementing Bayesian circular functional models (wrapped multivariate functional models)
- `Experiment.R`: Script to run experiments in the paper
- `Experiment-Summary.R`: Script to summarize experimental results to reproduce tables/figures in the paper 
- `PLV-original.RData`: Phase locking values (PLV) based on the original dataset obtained from [real EEG dataset available from the University of Cambridge Data Repository](https://www.repository.cam.ac.uk/handle/1810/252736}).
- `Data-Gaussian-noise.zip`: Zip file to include data with added Gaussian noise
- `Data-Gaussian-noise.zip`: Zip file to include data with added uniform noise


