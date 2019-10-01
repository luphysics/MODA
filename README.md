# MODA

## Introduction

MODA (Multiscale Oscillatory Dynamics Analysis) is a numerical toolbox developed by the
[Nonlinear & Biomedical Physics group](https://www.lancaster.ac.uk/physics/research/experimental-condensed-matter/nonlinear-and-biomedical-physics/) at [Lancaster University](https://www.lancaster.ac.uk/physics/) for analysing real-life time-series
that are assumed to be the output of some a priori unknown non-autonomous dynamical system,
and deriving important properties about this dynamical system from the time-series. It includes
methods both for analysing the recordings of a single signal over time, and for analysing a set
of recordings of multiple different signals over time. In particular, it has tools for analysing
bivariate time-series consisting of the simultaneous recordings of two different signals over time,
with a view to examining possible connections between the two signals.

The foundation of most of the methods in MODA is time-frequency analysis, which describes
the time-evolving frequency content of a signal. On the basis of this, one can examine whether
there is either mutual interaction between or common influence on a pair of signals, by way of
a “common” oscillatory component of possibly time-varying frequency, where “commonness” is
measured by coherence of phases. Similarly, one can examine whether there are interactions
between oscillatory components either within one signal or between two different signals,
using bispectral analysis based on the coherence of the sum of the phases with the extracted
phase from the harmonic frequency corresponding to the sum of the original frequencies; such
analysis can detect both nonlinear interactions between oscillatory components and also linear
interactions between nonlinear oscillatory components. Moreover, one can extract from a signal
and reconstruct oscillatory components with time-varying frequency, together with the evolution
of the “instantaneous phase” of such a component. From here, one can use dynamical Bayesian
inference to reconstruct approximations of a non-autonomous stochastic differential equation
describing the joint evolution of a pair of oscillatory components (of either the same signal or two different signals): under the approximation that phases of interacting oscillators evolve
according to the dynamics of coupled phase oscillators, one obtains over a series of time-windows
the “most likely” coupling function from within the span of a given number of orthonormal basis
functions on the 2-torus.


## Getting Started

### Preparation 

MODA requires a MATLAB installation. The latest version (R2019a) is recommended but not required.

The following MATLAB toolboxes should be installed:
- Signal Processing Toolbox                
- Statistics and Machine Learning Toolbox  
- Wavelet Toolbox         

You can check which toolboxes are currently installed by running the "ver" command in the MATLAB Command Window.

### Running MODA

To use MODA, download the code and place it in a desired location. In your file explorer, double-click "MODA.m" inside the MODA folder to open it with MATLAB. 

MODA can then be started using the "Run" button in the MATLAB editor.
