<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Example analysis for two interacting oscillators](#example-analysis-for-two-interacting-oscillators)
  - [Ridge extraction](#ridge-extraction)
    - [ECG](#ecg)
    - [Respiration](#respiration)
    - [Finding the frequencies](#finding-the-frequencies)
  - [Wavelet phase coherence](#wavelet-phase-coherence)
    - [Instantaneous heart rate and respiration](#instantaneous-heart-rate-and-respiration)
    - [Instantaneous heart rate and instantaneous respiration](#instantaneous-heart-rate-and-instantaneous-respiration)
  - [Dynamical Bayesian inference](#dynamical-bayesian-inference)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Example analysis for two interacting oscillators 

This analysis is based on two real time-series: 

- An ECG signal, showing heart rate. 
- A respiratory signal measured with a respiration belt. 

The signals were measured simultaneously and the recordings lasted approximately 30 minutes. The heart and the lungs can be considered as two oscillators, and we will investigate their interaction by using wavelet based time-frequency analysis. 

Both signals were resampled to a sampling frequency of **300Hz**, and they are saved in a **row-wise** configuration.

---

![Screenshot of the ECG signal](/docs/images/example_analysis/ecg_signal.png)

*A zoomed-in section of the ECG signal.*

---

![Screenshot of the respiration signal](/docs/images/example_analysis/resp_signal.png)

*A zoomed-in section of the respiration signal.*

---

## Ridge extraction

The first step of the analysis is to identify the frequency interval at which the oscillator is operating. This can be determined using ridge extraction. 

### ECG

- Start MODA. 
- In the launcher window, click "Ridge extraction and filtering". 
- In the toolbar, click "File" -> "Load time series". 
- Load the signal (`ECG.mat` from the `example_sigs/example_analysis` directory). 
- Enter the sampling frequency, 300Hz, when prompted. 
- Select row-wise as the orientation.
- Click "Calculate Transform".

---

![Screenshot of the wavelet transform of the ECG signal](/docs/images/example_analysis/ecg_wt.png)

*The wavelet transform of the ECG signal.*

---

The result can be used to identify a proper **minimum frequency**, **maximum frequency** and **frequency resolution**. These values can be entered in the "Transform Options" section in the upper right corner.

Higher values of frequency resolution (fr > 2) give good frequency resolution but poor time resolution; we found fr = 3 to be optimal for this signal. 

Since the heartrate is around 1Hz, the maximum frequency was set to 5Hz. Similarly, the minimum frequency was set to 0.095Hz.

Now the wavelet tranform can be calculated again using "Calculate Transform".

--- 

![Screenshot of the wavelet tranform of the ECG after adjusting the parameters.](/docs/images/example_analysis/ecg_wt2.png)

*The wavelet transform of the ECG signal after setting the minimum and maximum frequencies, and the frequency resolution.*

---

The region which should be analysed can be seen in the screenshot. It is a region with higher amplitudes, where the frequency is around an expected range (for heartrate, we expect a range of around 1Hz).

Mark a band around this region:

- Click "Mark region".
- Click on the plot, slightly below and above the desired region to set the frequencies. 

---

![Screenshot of the marked region in ridge extraction using the ECG signal](/docs/images/example_analysis/ecgmarkedregion.png)

*The marked region in ridge extraction.*

---

Now that the region has been marked:

- Click "Add marked region" to add it to the list of selected regions.
- Press "Extract ridges" in the bottom right corner.

The ridge will be extracted:

---

![Screenshot of the extracted ridge of the ECG signal](/docs/images/example_analysis/ecg_ridge.png)

*The extracted ridge of the ECG signal.*

---

![Screenshot of the result of ridge extraction on the ECG signal](/docs/images/example_analysis/ridge_ecg_zoom.png)

*The result of ridge extraction on the ECG signal (zoomed-in).*

---

The results can be saved as a mat file using "Save" -> "Save as .mat" in the toolbar. The ridge is the `Ridge_frequency` value in the saved matlab structure.

---

![Shows the structure saved from MODA.](/docs/images/Structure.png)

*The ridge data which was saved from MODA.*

---

### Respiration

The same process can be followed for the respiration signal. 

The wavelet transform is calculated:

---

![Screenshot of the wavelet transform of the respiration signal](/docs/images/example_analysis/resp_wt2.png)

*The wavelet transform of the respiration signal.*

---

For the respiration signal, a frequency resolution of fr = 1 is suitable. The minimum and maximum frequencies can also be left to their defaults.

The desired region is marked: 

---

![Screenshot of the band marking for the respiration signal](/docs/images/example_analysis/respmarkedregion.png)

*The marked region of the respiration signal.*

--- 

The ridge is extracted:

---

![Screenshot of the exracted ridge for the respiration signal](/docs/images/example_analysis/resp_ridge.png)

*The extracted ridge of the respiration signal.*

--- 

![Zoomed screenshot of the extracted ridge for the respiration signal](/docs/images/example_analysis/resp_ridge_zoomed.png)

*Zoomed-in screenshot of the extracted ridge for the respiration signal.*

---

As with the ECG signal, the result can be saved as a `.mat` file.

### Finding the frequencies

When performing ridge extraction, the frequency interval for each ridge can be obtained by using the data tips in "tools" when the figure is open, or by using the `max` and `min` functions on the saved `Ridge_frequency` arrays in MATLAB.

For example, in MATLAB the maximum and minimum frequencies can be obtained with the following code for ridge data saved as `ecg_ridge.mat`:

```matlab
data = load("ecg_ridge.mat")
f = data.Filtered_data

% Note: indices {1,1} are appropriate since only one signal pair was analysed.
ridge = f.Ridge_frequency{1,1}

fmin = min(ridge)
fmax = max(ridge)
```

The frequency intervals obtained for the ECG and respiration signals were:

| Signal | Min freq (Hz) | Max freq (Hz) | 
| --- | --- | --- |
| ECG | 0.91 | 1.3 | 
| Resp | 0.12 | 0.37 | 

## Wavelet phase coherence

Wavelet phase coherence can be used to investigate the similarities between two signals at particular frequencies. 

The wavelet phase coherence window loads a `.mat` file containing a single array containing two or more signals. This necessitates using MATLAB or another tool to create a data file which will be loaded.

In this section, we will analyse coherence between:

- Instantaneous heartrate and respiration.
- Instantaneous heartrate and instantaneous respiration.

The data files loaded in the Wavelet Phase Coherence window were produced using MATLAB; they are saved in the `example_sigs/example_analysis` folder as `ihrResp.mat` and `ihrirr.mat` respectively.

> **Note:** *Instantaneous* heartrate or respiration is the *extracted ridge* of the ECG or respiration signal respectively.

The frequency range you choose for this analysis should include the frequency intervals of both oscillators. Hence, an appropiate choice could be **0.12Hz to 1.3Hz**.

Coherence should always be tested for significance, and therefore we use surrogates. In this example we use **30 Fourier Transform (FT) surrogates**.

> **Note:** For more information about surrogates, please see the [User Manual](/User%20Manual.pdf)).

### Instantaneous heart rate and respiration

- Open MODA.
- Click "Wavelet Phase Coherence".
- Click "File" -> "Load time series".
- Browse to, and select, the first data file (`example_sigs/example_analysis/ihrResp.mat`).
- Enter the sampling frequency (300Hz).
- Enter the orientation (row-wise).
- Enter the minimum and maximum frequencies.
- Enter `30` in the "Surrogate Count" box.
- Choose `FT` in the surrogate "Method" dropdown.
- Click "Coherence - all".

> **Note:** Due to the number of surrogates, this calculation may take some time.

The wavelet phase coherence is calculated:

---

![Figure showing the Wavelet Phase Coherence in MODA, ECG ridge and respiration signal.](/docs/images/example_analysis/ihrrespcoherenceav.png)

*Screenshot showing the wavelet phase coherence between the ECG ridge and respiration signal.*

---

### Instantaneous heart rate and instantaneous respiration

The above process can be repeated, using the `ihrirr.mat` data file.

---

![Figure showing the Wavelet Phase Coherence in MODA, ECG ridge and respiration ridge.](/docs/images/example_analysis/ihrirrcoherenceav.png)

*Screenshot showing the wavelet phase coherence between the ECG ridge and respiration ridge.*

---

## Dynamical Bayesian inference

Now press "Dynamical Bayesian inference" in the MODA GUI, and load the matrix with two signals. 

The frequency ranges are now specified for each oscillator. In the box "Freq range 1" write the frequency interval for the first signal (0.8 to 1.3Hz for ECG), and in freq range 2 write the frequency interval for the second signal (0.13 to 0.29Hz for respiration). Then choose a number of surrogates (we have used 50 in this example), and press "Add parameter set". After this you can press "Calculate". 

---

![Figure showing the dynamical Bayesian inference in MODA, ECG ridge and respiration signal.](/docs/images/BayesianIHRResp.png)

*Screenshot showing dynamical Bayesian inference in MODA, using the ECG ridge and respiration signal.*

---

![Figure showing the dynamical Bayesian inference in MODA, ECG ridge and respiration ridge.](/docs/images/BayesianIHRIRR.png)

*Screenshot showing dynamical Bayesian inference in MODA, using the ECG ridge and respiration ridge.*

---

The coupling function can be viewed by clicking the dropdown in the bottom-left corner, and selecting "Coupling function".

---

![Screenshot showing the coupling function with the ECG ridge and respiration signal.](/docs/images/CF_IHRIRR.png)

*Screenshot showing the coupling function with the ECG ridge and respiration signal.*

---

![Screenshot showing the coupling function with the ECG ridge and respiration ridge.](/docs/images/CF_IHRResp.png)

*Screenshot showing the coupling function with the ECG ridge and respiration ridge.*

---
