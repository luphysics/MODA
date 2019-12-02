## User Guide

This guide describes how to get started with MODA. More in-depth instructions about using MODA's features can be found in the user manual, `User Manual.pdf`, which will be downloaded with the source code. 

### Preparation 

MODA requires a MATLAB installation. MATLAB R2017a or higher is required, but newer versions are recommended.

The following MATLAB toolboxes are needed:
- Signal Processing Toolbox                
- Statistics and Machine Learning Toolbox  
- Wavelet Toolbox         

You can check which toolboxes are currently installed by running the `ver` command in the MATLAB Command Window.

### Downloading MODA

There are two methods to download the code: as a zip file, or by cloning the repository with Git. The advantage of cloning with Git is that you can easily update MODA in-place by running `git pull` in the terminal, preserving any added files such as shortcuts, instead of downloading a new zip file.

If you prefer the zip method:

- [Click here](https://github.com/luphysics/MODA/zipball/master) to download the .zip file. 
- Extract the zip file to a desired location.
- For simplicity of instructions, rename the folder to `MODA`. 

If you prefer the Git method:

- [Install Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
- Open a terminal in a desired folder and run `git clone https://github.com/luphysics/MODA.git`.
- The code will download as a folder named `MODA`.
- Whenever you want to update, open a terminal in `MODA` and run `git pull`.

### Running MODA

In your file explorer, double-click `MODA.m` inside the `MODA` folder to open it with MATLAB. After the MATLAB window opens, press F5 or click the "Run" button to start MODA.

If the following dialog appears, click "Change Folder":

![Screenshot of the dialog stating that MODA.m is not in the current MATLAB path.](/docs/images/change_folder.png)

> **Note:** You may need to click inside the section displaying the contents of `MODA.m` for the "Run" button to appear.

### Importing data

When MODA opens, try clicking "Time-Frequency Analysis". After the window opens, go to "File" -> "Load time series" in the top left of the window. Using the file browser dialog, select a `.csv` or `.mat` file.

> **Tip:** There are some example signals in the `example_sigs` folder. Try `example_sigs/6signals_10Hz.mat`.

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

![Screenshot of the time-frequency window after data is loaded.](/docs/images/timefrequency_empty.png)