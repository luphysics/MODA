<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [User Guide](#user-guide)
  - [Getting started](#getting-started)
    - [Requirements](#requirements)
    - [Downloading MODA](#downloading-moda)
    - [Running MODA](#running-moda)
    - [Importing time-series](#importing-time-series)
    - [Example - Time-Frequency Analysis](#example---time-frequency-analysis)
    - [Truncating signals](#truncating-signals)
    - [Plotting and saving](#plotting-and-saving)
    - [Large arrays and downsampling](#large-arrays-and-downsampling)
  - [Launcher window](#launcher-window)
  - [Time-Frequency Analysis](#time-frequency-analysis)
  - [Ridge Extraction and Filtering](#ridge-extraction-and-filtering)
    - [Negative frequencies](#negative-frequencies)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# User Guide

This guide is aimed at users wishing to set up and use MODA. The [User Manual](/User%20Manual.pdf) provides a more in-depth explanation of MODA's functionality.

If you're interested in modifying or contributing to the program, you may find the [Developer Guide](/docs/developer-guide.md) useful.

## Getting started

### Requirements

MATLAB R2017a or higher is required, but newer versions are recommended.

The following MATLAB toolboxes are needed:
- Signal Processing Toolbox                
- Statistics and Machine Learning Toolbox  
- Wavelet Toolbox         

You can check which toolboxes are currently installed by running the `ver` command in the MATLAB Command Window.

### Downloading MODA

- [Click here](https://github.com/luphysics/MODA/zipball/master) to download the code as a .zip file. 
- Extract the zip file to a desired location.
- For simplicity of instructions, rename the folder to `MODA`. 

### Running MODA

In your file explorer, double-click `MODA.m` inside the `MODA` folder to open it with MATLAB. After the MATLAB window opens, press `F5` or click the "Run" button to start MODA.

> **Note:** You may need to click inside the section displaying the contents of `MODA.m` for the "Run" button to appear.

If the following dialog appears, click "Change Folder":

![Screenshot of the dialog stating that MODA.m is not in the current MATLAB path.](/docs/images/change_folder.png)

The launcher window will then open:

![Screenshot of MODA's launcher window.](/docs/images/launcher_window.png)

### Importing time-series

In MODA, a time-series is a series of recorded values, where the sampling frequency - the frequency at which the recordings were made - is known.

MODA can analyse multiple signals, provided that all signals have the same duration and sampling frequency.

To import time-series into MODA, they must be saved in a compatible format: 

- The file type must be a `.mat` file or `.csv` file. 
- The file must contain a **single array, whose entries are all a single real number**. Each row or column of the array corresponds to a different time-series. 
- For windows which inspect pairs of signals (for example, dynamical Bayesian inference), the number of time-series should be even. *If there are an odd number of time-series, the last one will be removed and pairs will be formed from the remaining time-series.*
- For wavelet bispectrum analysis, there must be only two time-series; this is because the bispectrum computations take a prohibitively long time.

> **Note:** The sampling frequency is entered in the user interface, rather than specified in the file.

**If the array loaded into MODA is extremely large, it may run slowly or crash.** A potential way to overcome this is to downsample the signals.

### Example - Time-Frequency Analysis

When MODA opens, try clicking "Time-Frequency Analysis". After the window opens, go to `File` -> `Load time series` in the top left of the window. Using the file browser dialog, select a `.csv` or `.mat` file.

> **Tip:** There are some example signals in the `example_sigs` folder. Try `example_sigs/6signals_10Hz.mat` (a row-wise signal).

After the file is selected, a dialog will appear: 

![Screenshot of the sampling frequency dialog.](/docs/images/sampling_frequency.png)

After entering the sampling frequency in Hz (for example, "10"), another dialog will appear:

![Screenshot of the data orientation dialog.](/docs/images/data_orientation.png)

This dialog asks whether the orientation of the data is row-wise or column-wise.

---

###### Row-wise data

With row-wise data, **each row contains one signal**.

```
| Signal 1, Value 1 | Signal 1, Value 2 | Signal 1, Value 3 | 
| Signal 2, Value 1 | Signal 2, Value 2 | Signal 2, Value 3 | 
| Signal 3, Value 1 | Signal 3, Value 2 | Signal 3, Value 3 | 
| Signal 4, Value 1 | Signal 4, Value 2 | Signal 4, Value 3 | 
```

###### Column-wise data

With column-wise data, **each column contains one signal**.

```
| Signal 1, Value 1 | Signal 2, Value 1 | Signal 3, Value 1 | 
| Signal 1, Value 2 | Signal 2, Value 2 | Signal 3, Value 2 | 
| Signal 1, Value 3 | Signal 2, Value 3 | Signal 3, Value 3 | 
| Signal 1, Value 4 | Signal 2, Value 4 | Signal 3, Value 4 | 
```

---

After selecting the orientation, the data will be loaded and the first signal will be plotted at the top of the window.

> **Note:** If the wrong orientation is selected, MODA may freeze because it interprets the data as a very large number of short signals.

![Screenshot of the time-frequency window after data is loaded.](/docs/images/timefrequency_empty.png)

The imported signals are listed in the "Select data" section which, in the screenshot, is outlined in red. If the window's functionality demands that signals are processed as pairs, then each item in the list will represent a particular signal pair.

When a signal or signal pair is selected, it (and its associated results, if any have been calculated) will be plotted. 

### Truncating signals

You may wish to **analyse only a portion of the recorded signal(s).** When the signals have been imported, you can do this by zoom-in magnifying glass and then using a click-and-drag action on the signal plot, to zoom to a rectangular region.

![Screenshot demonstrating the click-and-drag action to zoom to a region.](/docs/images/zoom_signal.png)

> **Tip:** When selecting a portion of the signal, the signal will zoom to the horizontal boundaries of the selection; the vertical zoom will be reset.

Alternatively, it is possible to set the by entering values in the "Xlim" field and clicking the "Refresh" button.

![Screenshot demonstrating the "Xlim" field.](/docs/images/xlim.png)

After the signal has been truncated, the plot of the main results will not change until the calculations are repeated. 

> **Note:** For results involving frequency or time-frequency analysis, **the minumum frequency for which results can be displayed** increases as the signal is truncated.

### Plotting and saving

**To save one of the graphs currently displayed** in the window, click `Plot` in the toolbar (to the right of `File`) and then select the appropriate option; this will open the plot in a new window as a MATLAB figure.

This has some benefits:

- The figure can be saved in a variety of formats.
- All of MATLAB's features for analysing and manipulating a figure are available.
- It will not be overwritten with the results of a new calculation.

When saving figures, the following formats are recommended:

- For a plot displaying a function of a single variable, a vector image format such as `.svg`, `.pdf` or `.eps`.
- For a plot displaying a function of two variables, a scalar image format such as `.png`.

#### Saving the current session

To save the current session, click `Save`->`Save session` in the toolbar. This will save the current session as a `.mat` file, allowing you to return to it later by using `File`->`Load previous session`.

#### Saving data

Many of the actual values displayed in the plots can be saved as files by clicking `Save`->`Save as .mat` or `Save`->`Save as .csv` in the toolbar.

`.mat` files can be loaded and interpreted in MATLAB or Python; `.csv` files are useful for Microsoft Excel.

> :warning: When saving as a .csv file, avoid changing the file extension in the filename, or changing the file type in the "Save as type" dropdown.

Some numerical values will be saved as NaN ("not a number"), which represents that a value has not been computed. In time-frequency analysis, meaningful results cannot be found very close to the start/end of the signal; if "cut edges" is enabled, these values will be saved as NaN.

> :warning: When saving data from Ridge Extraction & Filtering for the purpose of loading into Dynamical Bayesian Inference, results must be saved as a `.mat` file.

> **Note:** Files created by MODA are often very large; files may not save immediately, and opening a partially saved file could corrupt it.

### Large arrays and downsampling

If an extremely large array is loaded into MODA, it could run very slowly or crash.

Arrays may be overly large for the following reasons:

1. The number of signals stored is large.
2. The signals were measured over a very long time.
3. The signals were measured with a very high sampling frequency.

Issue (1) can be solved by splitting up the array into multiple files, each containing fewer signals, and analysing each file individually.

Issue (2) can be solved by splitting up the time interval over which the signals were measured into smaller intervals (preferably with a small overlap, since time-frequency analysis cannot be performed too close to the start/end of a signal), saving as smaller arrays, and analysing each array.

If the sampling frequency is many times higher (e.g. ~10 times higher) than the largest potential frequency of interest within the signals, issue (3) can be addressed by **downsampling**.

## Launcher window

When MODA launches, it will open the launcher window:

![Screenshot of MODA's launcher window.](/docs/images/launcher_window.png)

The launcher window provides buttons to open the main windows, which are each described by in the sections below.

## Time-Frequency Analysis

The Time-Frequency Analysis window opens as follows:

![Screenshot of the time-frequency analysis window.](/docs/images/window_tf_blank.png)

After loading a data file, the time-series and pre-processed time-series will be plotted:

![Screenshot of the time-frequency analysis window after loading a data file.](/docs/images/window_tf_loaded.png)

The basic options are as follows:

- The dropdown titled "Transform" can be used to select the transform type as **wavelet transform** (WT) or **windowed Fourier transform** (WFT).
- The range of frequencies to analyse can be modified by entering values in the "Transform Options" section. If values are chosen, the "Max Freq" must be **at most half of the sampling frequency** of the signal. If left blank, the maximum frequency will default to half of the sampling frequency. **For the WFT, the minimum frequency must be specified**.
- The frequency resolution can also be entered in the "Transform Options" section. When left blank, the default value will be 1.

The "Advanced Options" section can be used as follows:
- The "WT/WFT Type" dropdown allows different wavelet types to be chosen for the WT, and different window types for the WFT.
- Pre-processing can be enabled/disabled. When pre-processing is enabled, the signal is "de-trended" by subtracting a best-fit cubic approximation. This may be useful for identifying oscillatory components of the signal, because non-oscillatory trends in the signal will affect the time-frequency representation. If a frequency interval has been specified, the frequencies outside the range will be filtered out. If pre-processing is disabled, the transform is simply applied to the original signal.
- Since time-frequency analysis methods rely on taking some “spread” about each moment in
time (as determined by the parameter fr ), one cannot meaningfully assign amplitudes and
phases to points in time-frequency space that are very close to the start of the signal or
very close to the end of the signal. When the option Cut Edges is on, MODA will respect
this fact by simply not attempting to calculate the selected transform at points where
this cannot sensibly be done; in the graphical representation of the transform that MODA
produces, those regions in the time-frequency space where the transform is not computed
will be left white. However, if the user decides to turn off the Cut Edges option, then
MODA will show the results as defined for the “zero-padding extension” of the original
signal, where the signal is simply taken to be equal to 0 before the start and after the end
of the actual recorded segment.

To perform a transform:

- For only the currently selected signal, click "Transform Single".
- For all signals, click "Transform All".

After performing a transform, the window displays the result:

![Screenshot of the time-frequency analysis window after calculating a transform.](/docs/images/window_tf_wt.png)

## Ridge Extraction and Filtering

### Negative frequencies

In some cases, negative frequencies may be present in the result of ridge extraction.

When this occurs, a dialog will be shown:

![Screenshot of the dialog which appears when negative frequencies are present in the result of ridge extraction.](/docs/images/dialog_negative_freq.png)

Ridge extraction uses one of two reconstruction algorithms: *ridge* or *direct* reconstruction. 

Direct reconstruction is the default. Direct reconstruction for frequency is inferred directly from the exact component reconstruction; for a signal such as `a(t)*cos(phi(t))` with slowly varying amplitude and frequency, it is capable of calculating the *exact* amplitude, phase and frequency, whereas the ridge reconstruction algorithm cannot. 

However, direct reconstruction is more susceptible to noise and interference between multiple components, which can result in negative frequencies appearing for certain signals.

Overall, the direct method is better for signals where the components are relatively clear. The ridge method is more robust for signals which especially suffer from noise and interference.

To address problems with negative frequencies, try one or more of the following:

- Switching to the 'Lognorm' wavelet.
- Using a narrower frequency interval for ridge extraction.
- Using a higher frequency resolution.

If none of these solve the issue, you could try switching to the ridge reconstruction algorithm. 

#### Switching to ridge reconstruction

Switching to ridge reconstruction involves changing a single line of code.

In the MATLAB Editor, open the file at `allguis/guis/filtering-master/Functions/MODAridge_filter.m`. 

Find this section (around line 80):

```matlab
% TO SWITCH FROM "direct" TO "ridge" RECONSTRUCTION, CHANGE 'direct' TO 'ridge' ON THE NEXT LINE.
[handles.bands_iamp{j,k},handles.bands_iphi{j,k},handles.bands_freq{j,k}] = rectfr(tfsupp,WT,freqarr,wopt,'direct');
```

Then change the `'direct'` parameter to `'ridge'` in the `rectfr` function call.