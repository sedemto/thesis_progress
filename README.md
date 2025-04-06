# PHAIN (PHase-aware Audio INpainter)
**Tomoro Tanaka (Department of Intermedia Art and Science, Waseda University, Tokyo, Japan)**\
[![DOI](https://zenodo.org/badge/690949058.svg)](https://zenodo.org/doi/10.5281/zenodo.10818129)

This README file describes the MATLAB codes provided to test, analyze, and evaluate our proposed method, PHAINs.\
PHAIN is an audio inpainting method introduced in the following paper
>[1] Tomoro Tanaka, Kohei Yatabe, and Yasuhiro Oikawa, "PHAIN: Audio Inpainting via Phase-aware Optimization with Instantaneous Frequency".

## Requirements
The codes were developed in MATLAB version R2023a and have been tested in R2023a and R2022b.\
Some functions rely on 

1. MathWorks Toolbox: You are kindly requested to download some of them, such as 'Signal Processing Toolbox'.

2. Toolbox available online: This is available online under the [MIT license](https://opensource.org/licenses/mit-license.php).

- Instantaneous phase correction of complex spectrogram\
  This folder contains a set of MATLAB codes written as a supplementary material of the following tutorial paper explaining phase-related property of spectrograms:
  I already installed it so you can easily execute the codes. Plaese refer to https://doi.org/10/c3qb.

  >[2] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa, "Representation of complex spectrogram via phase conversion," Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

## Usage
Execute `main.m` and run PHAINs (the parameters and data directory can be modified by yourself).

## Directory Structure

- `PHAIN`
  - `PHAINmain.m`.
  - `CP.m` is the Chambolle-Pock algorithm.

- `dataset` is the folder to contain test data, where any files can be downloaded.

- `utils`
  - `projGamma.m`.
  - `shortenForDGT.m`.
  - `calcSNR.m`.

- `phase_correction` is the toolbox mentioned above (see Requirements 2).

- `utils` contains some useful functions


## License
See the file named `LICENSE`.
