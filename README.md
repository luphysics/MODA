<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [MODA](#moda)
  - [Purpose](#purpose)
- [User Guide](#user-guide)
  - [Requirements](#requirements)
  - [Downloading MODA](#downloading-moda)
  - [Running MODA](#running-moda)
  - [Importing time-series](#importing-time-series)
- [References](#references)
  - [Example applications](#example-applications)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[![DOI](https://zenodo.org/badge/194114858.svg)](https://zenodo.org/badge/latestdoi/194114858)

# MODA

MODA (Multiscale Oscillatory Dynamics Analysis) is a numerical toolbox developed by the
[Nonlinear & Biomedical Physics group](https://www.lancaster.ac.uk/physics/research/experimental-condensed-matter/nonlinear-and-biomedical-physics/) at [Lancaster University](https://www.lancaster.ac.uk/physics/) and the Nonlinear Dynamics and Synergetic Group at the Faculty of Electrical
Engineering, University of Ljubljana, Slovenia under the supervision of Aneta Stefanovska.

To get started, please see the [User Guide](#user-guide).

> **Note:** A Python implementation of MODA, [PyMODA](https://github.com/luphysics/PyMODA), is currently in development. PyMODA does not require a MATLAB license.

## Purpose

MODA is designed for analysing real-life time-series
that are assumed to be the output of some *a priori* unknown non-autonomous dynamical system,
and deriving important properties about this dynamical system from the time-series. It includes
methods both for analysing the recordings of a single signal over time, and for analysing a set
of recordings of multiple different signals over time. In particular, it has tools for analysing
bivariate time-series consisting of the simultaneous recordings of two different signals over time,
with a view to examining possible connections between the two signals.

# User Guide

This guide is aimed at users wishing to set up and use MODA. The [User Manual](/User%20Manual.pdf) provides a more in-depth explanation of MODA's functionality.

If you're interested in modifying or contributing to the program, you may find the [Developer Guide](/docs/developer-guide.md) useful.

## Requirements

MATLAB R2017a or higher is required, but newer versions are recommended.

The following MATLAB toolboxes are needed:
- Signal Processing Toolbox                
- Statistics and Machine Learning Toolbox  
- Wavelet Toolbox         

You can check which toolboxes are currently installed by running the `ver` command in the MATLAB Command Window.

## Downloading MODA

- [Click here](https://github.com/luphysics/MODA/zipball/master) to download the code as a .zip file. 
- Extract the zip file to a desired location.
- For simplicity of instructions, rename the folder to `MODA`. 

## Running MODA

In your file explorer, double-click `MODA.m` inside the `MODA` folder to open it with MATLAB. After the MATLAB window opens, press `F5` or click the "Run" button to start MODA.

> **Note:** You may need to click inside the section displaying the contents of `MODA.m` for the "Run" button to appear.

If the following dialog appears, click "Change Folder":

![Screenshot of the dialog stating that MODA.m is not in the current MATLAB path.](/docs/images/change_folder.png)

The launcher window will then open:

![Screenshot of MODA's launcher window.](/docs/images/launcher_window.png)

## Importing time-series

In MODA, a time-series is a series of recorded values, where the sampling frequency - the frequency at which the recordings were made - is known.

MODA can analyse multiple signals, provided that all signals have the same duration and sampling frequency.

To import time-series into MODA, they must be saved in a compatible format: 

- The file type must be a `.mat` file or `.csv` file. 
- The file must contain a **single array, whose entries are all a single real number**. Each row or column of the array corresponds to a different time-series. 
- For windows which inspect pairs of signals (for example, dynamical Bayesian inference), the number of time-series should be even. *If there are an odd number of time-series, the last one will be removed and pairs will be formed from the remaining time-series.*
- For wavelet bispectrum analysis, there must be only two time-series; this is because the bispectrum computations take a prohibitively long time.

> **Note:** The file should only contain the values of the time-series, because the sampling frequency is entered in the user interface.

**If the array loaded into MODA is extremely large, it may run slowly or crash.** A potential way to overcome this is to downsample the signals.

### Example

When MODA opens, try clicking "Time-Frequency Analysis". After the window opens, go to `File` -> `Load time series` in the top left of the window. Using the file browser dialog, select a `.csv` or `.mat` file.

> **Tip:** There are some example signals in the `example_sigs` folder. Try `example_sigs/6signals_10Hz.mat` (a row-wise signal).

After the file is selected, a dialog will appear: 

![Screenshot of the sampling frequency dialog.](/docs/images/sampling_frequency.png)

After entering the sampling frequency in Hz (for example, "10"), another dialog will appear:

![Screenshot of the data orientation dialog.](/docs/images/data_orientation.png)

This dialog asks whether the data is row-wise or column-wise.

---

#### Row-wise data

With row-wise data, each row corresponds to a different signal.

```
| Signal 1, Value 1 | Signal 1, Value 2 | Signal 1, Value 3 | 
| Signal 2, Value 1 | Signal 2, Value 2 | Signal 2, Value 3 | 
| Signal 3, Value 1 | Signal 3, Value 2 | Signal 3, Value 3 | 
| Signal 4, Value 1 | Signal 4, Value 2 | Signal 4, Value 3 | 
```

#### Column-wise data

With column-wise data, each column corresponds to a different signal.

```
| Signal 1, Value 1 | Signal 2, Value 1 | Signal 3, Value 1 | 
| Signal 1, Value 2 | Signal 2, Value 2 | Signal 3, Value 2 | 
| Signal 1, Value 3 | Signal 2, Value 3 | Signal 3, Value 3 | 
| Signal 1, Value 4 | Signal 2, Value 4 | Signal 3, Value 4 | 
```

---

After selecting the orientation, the data will be loaded and the first signal will be plotted at the top of the window.

> **Note:** If the wrong orientation is selected, MODA may freeze while attempting to load the data as a large number of small signals.

![Screenshot of the time-frequency window after data is loaded.](/docs/images/timefrequency_empty.png)

The signals are listed in the "Select data" section (outlined in red), and each signal will be plotted when selected.

# References

#### Overview
1. J Newman, G Lancaster and A Stefanovska, “Multiscale Oscillatory Dynamics
Analysis”, v1.01, User Manual, 2018.
2. P Clemson, G Lancaster, A Stefanovska, “Reconstructing time-dependent dynamics”, *Proc IEEE*
**104**, 223–241 (2016).
3. P Clemson, A Stefanovska, “Discerning non-autonomous dynamics”, *Phys Rep* **542**, 297-368 (2014).

#### Time-Frequency Analysis
1. D Iatsenko, P V E McClintock, A Stefanovska, “Linear and synchrosqueezed time-frequency
representations revisited: Overview, standards of use, resolution, reconstruction, concentration, and
algorithms”, *Dig Sig Proc* **42**, 1–26 (2015).
2. P Clemson, G Lancaster, A Stefanovska, “Reconstructing time-dependent dynamics”, *Proc IEEE*
**104**, 223–241 (2016).
3. G Lancaster, D Iatsenko, A Pidde, V Ticcinelli, A Stefanovska, “Surrogate data for hypothesis testing of
physical systems”, *Phys Rep* **748**, 1–60 (2018).

#### Wavelet Phase Coherence
1. Bandrivskyy A, Bernjak A, McClintock P V E, Stefanovska A, “Wavelet phase coherence analysis:
Application to skin temperature and blood flow”, *Cardiovasc Engin* **4**, 89–93 (2004).
2. Sheppard L W, Stefanovska A, McClintock P V E, “Testing for time-localised coherence in bivariate
data”, *Phys. Rev. E* **85**, 046205 (2012).

#### Ridge Extraction & Filtering
1.  D Iatsenko, P V E McClintock, A Stefanovska, “Nonlinear mode decomposition: A noise-robust,
adaptive decomposition method”, *Phys Rev E* **92**, 032916 (2015).
2. D Iatsenko, P V E McClintock, A Stefanovska, “Extraction of instantaneous frequencies from ridges in
time-frequency representations of signals”, *Sig Process* **125**, 290–303 (2016).

#### Wavelet Bispectrum Analysis
1. J Jamšek, A Stefanovska, P V E McClintock, “Wavelet bispectral analysis for the study of interactions
among oscillators whose basic frequencies are significantly time variable”, *Phys Rev E* **76**, 046221
(2007).
2. J Jamšek, M Paluš, A Stefanovska, “Detecting couplings between interacting oscillators with
time-varying basic frequencies: Instantaneous wavelet bispectrum and information theoretic approach”,
*Phys Rev E* **81**, 036207 (2010).
3. J Newman, A Pidde, A Stefanovska, “Defining the wavelet bispectrum”, submitted (2019).

#### Dynamical Bayesian Inference
1. V N Smelyanskiy, D G Luchinsky, A Stefanovska, P V E McClintock, “Inference of a nonlinear stochastic model of the cardiorespiratory
interaction”, *Phys Rev Lett* **94**, 098101 (2005).
2. T Stankovski, A Duggento, P V E McClintock, A Stefanovska, “Inference of time-evolving coupled dynamical systems in the presence of noise”,
*Phys Rev Lett* **109**, 024101 (2012).
3. T Stankovski, A Duggento, P V E McClintock, A Stefanovska, “A tutorial on time-evolving dynamical Bayesian inference”, *Eur Phys J – Special
Topics* **223**, 2685-2703 (2014).
4. T Stankovski, T Pereira, P V E McClintock, A Stefanovska, “Coupling functions: Universal insights into dynamical interaction mechanisms”, *Rev
Mod Phys* **89**, 045001 (2017).
5. Special issue of the *Philos Trans Royal Soc A* (2019) with contributions by Kuramoto and others.

## Example applications

#### Wavelet Phase Coherence
1. Sheppard L W, Vuksanović V, McClintock P V E, Stefanovska A, Oscillatory dynamics of
vasoconstriction and vasodilation identified by time-localized phase coherence *Phys Med Biol*
**56**, 3583–3601 (2011).
2. A Bernjak, J Cui, S Iwase, T Mano, A Stefanovska, D L Eckberg, “Human sympathetic outflows to skin
and muscle target organs fluctuate concordantly over a wide range of time-varying frequencies”, *J
Physiol* **590**, 363–375 (2012).
3. P Kvandal, L Sheppard, S A Landsverk, A Stefanovska, K A Kirkebøen, “Impaired cerebrovascular
reactivity after acute traumatic brain injury can be detected by wavelet phase coherence analysis of the
intracranial and arterial blood pressure signals”, *J Clin Monit Comput* **27**, 375-383 (2013).

#### Ridge Extraction & Filtering
1. D Iatsenko, A Bernjak, T Stankovski, Y Shiogai, P J Owen-Lynch, P B M Clarkson, P V E McClintock,
A Stefanovska, “Evolution of cardiorespiratory interactions with age”, *Phil Trans R Soc A* **371**,
20110622 (2013).
2. V Ticcinelli, T Stankovski, D Iatsenko, A Bernjak, A E Bradbury, A R Gallagher, P B M Clarkson, P V
E McClintock, A Stefanovska, “Coherence and coupling functions reveal microvascular impairment in
treated hypertension”, *Front Physiol* **8**, 749 (2017).
3. YA Abdulhameed, G Lancaster, PVE McClintock, A Stefanovska, “On the suitability of laser-Doppler
flowmetry for capturing microvascular blood flow dynamics from darkly pigmented skin”, *Physiol Meas*,
**40**, 074005 (2019).

#### Wavelet Bispectrum Analysis
1. J Jamšek, A Stefanovska, P V E McClintock, “Nonlinear cardio-respiratory interactions revealed by
time-phase bispectral analysis”, *Phys Medicine Biol* **49**, 4407 (2004).

#### Dynamical Bayesian Inference
1. B Musizza, A Stefanovska, P V E McClintock, M Paluš, J Petrovčič, S Ribarič, F F Bajrović, “Interactions between cardiac, respiratory and
EEG-delta oscillations in rats during anaesthesia”, *J Physiol* **580** 315–326 (2007).
2. T Stankovski, V Ticcinelli, P V E McClintock, A Stefanovska, “Coupling functions in networks of oscillators”, *New J Phys* **17**, 035002 (2015).
3. T Stankovski, S Petkoski, J Ræder, A F Smith, P V E McClintock, A Stefanovska, “Alterations in the coupling functions between cortical and
cardio-respiratory oscillations due to anaesthesia with propofol and sevoflurane”, *Philos Trans Royal Soc A* **374**, 20150186 (2016).
4. V Ticcinelli, T Stankovski, D Iatsenko, A Bernjak, A E Bradbury, A R Gallagher, P B M Clarkson, P V E McClintock, A Stefanovska, “Coherence
and coupling functions reveal microvascular impairment in treated hypertension”, *Front Physiol* **8**, 749 (2017).
