# mosickosci – Posturography, motion sickness susceptibility and multisensory re-weighting 

## Introduction

The present repository contains the data and the code for analysing the
data related to an experiment involving the relationship between
posturographic measurements and motion sickness susceptibility. This
relationship is investigated through the analysis of the spectrum of the
center of pressure signal when the participant keeps balance on a moving
platform, with eyes either open or closed.

For more details, see the [OSF
page](https://osf.io/buk74/?view_only=a95b84f14cfe4669a3f1fa1334def588) for
the project.


## Installation

The code can be obtained with the following command:

```
git clone https://github.com/Posture-Motion-Sickness/mosickosci
```

N.B.: Developers with write access right to the repository should do,
instead:

```
git clone git@github.com:Posture-Motion-Sickness/mosickosci.git mosickosci
```

## Running the analysis and producing the figures

The scripts are written in [R](https://www.r-project.org/). The needed R
packages are automatically installed as required by the scripts. For
generating some PDF files with te composite results, the system program
`make` is also required.


## Description of the contents of this repository

### Data

* [data/mssq.csv](data/mssq.csv): Participant data from the Motion Sickness
  Susceptibility Questionnaire 
* [data/session-order.csv](data/session-order.csv): Order of trails for
  each participant
* [data/sNNN-YY-FF.brute](data/): Individual, per-condition posturographic
  data, where `NNN` is the participant number, `YY` is the vision condition
  (`yo`: open eyes, `yf`: closed eyes), and `FF` is the platform
  oscillation frequency (`f1`: 0.1 Hz, `f2`: 0.2 Hz, and `f4`: 0.4 Hz)


### Scripts for the analysis

* [config.r](analysis/config.r): General configuration and utility functions
* [raw-to-pca.r.r](analysis/raw-to-pca.r): Compute the principal
  component analysis on the PSD from the raw data 
* [pc1-variance.r](analysis/pc1-variance.r): Descriptive analysis of the PC1 score
* [pc1-score-analysis.r](analysis/pc1-score-analysis.r): Analysis of the
  temporal evolution of the PC1 score
* [spectral-balance.r](analysis/spectral-balance.r): Spectral balance
  analysis of the PC1 eigenvectors


## Licensing conditions

The files in this repository are made available under the conditions of the
[GNU Public License, version 3 or later](COPYING). No warranties. The study
related to this repository has been submitted for publication. If you use
the software or the data of this repository, please give credit. Stay tuned
for more information on the eventual acceptation for publication of the
manuscript .


## Authors

* Copyright © 2021 [Rafael Laboissière](https://github.com/rlaboiss)
* Copyright © 2021 [Merrick Dida](https://github.com/DidaMer)
