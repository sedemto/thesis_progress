%% Display explanation of the codes
disp ' '
disp ' '
disp '   Demonstration of instantaneous phase correction of complex spectrogram.'
disp ' '
disp ' '
disp '1. INTRODUCTION'
disp ' '
disp 'This script serves as a README document and run script in Code Ocean. If you'
disp 'have your own MATLAB, it should be easier to run each demo code separately'
disp 'after downloading the folder named "DGT+iPC".'
disp ' '
disp 'The set of MATLAB codes was written as a supplementary material of the following'
disp 'paper explaining phase-related property of the discrete Gabor transform (DGT):'
disp ' '
disp ' [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,'
disp '     "Representation of complex spectrogram via phase conversion,"'
disp '     Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)'
disp ' '
disp 'These codes provide some simple methods to perform the following things:'
disp ' '
disp ' (1) Calculating spectrogram of a real-valued signal;'
disp ' (2) Designing a synthesis window for perfectly reconstructing the signal;'
disp ' (3) Calculating the inverse transform reconstructing the signal from its spectrogram;'
disp ' (4) Calculating instantaneous frequency at each time-frequency bin;'
disp ' (5) Obtaining derivative of a window for computing instantaneous frequency; and'
disp ' (6) Calculating instantaneous-phase-corrected (iPC) spectrogram explained in [1].'
disp ' '
disp ' '
disp '2. INCLUDED FUNCTIONS'
disp ' '
disp 'The folder "DGT+iPC" contains the following MATLAB functions:'
disp ' '
disp '  Executable demo files:'
disp '   - demo1_DGTusage.m'
disp '   - demo2_windowUsage.m'
disp '   - demo3_IFandIPC.m'
disp ' '
disp '  Window function utilities:'
disp '   - generalizedCosWin.m'
disp '   - calcCanonicalDualWindow.m'
disp '   - calcCanonicalTightWindow.m'
disp '   - numericalDiffWin.m'
disp ' '
disp '  Discrete Gabor transform:'
disp '   - DGT.m'
disp '   - invDGT.m'
disp '   - zeroPaddingForDGT.m'
disp ' '
disp '  Instantaneous phase correction:'
disp '   - calcInstFreq.m'
disp '   - instPhaseCorrection.m'
disp '   - invInstPhaseCorrection.m'
disp ' '
disp 'In addition to these files, a faster implementation of DGT and its inverse is'
disp 'also included in the subfolder for practical use:'
disp ' '
disp '  Executable demo file:'
disp '   - demo_FDGTusage.m'
disp ' '
disp '  Faster DGT implementation:'
disp '   - FDGT.m'
disp '   - invFDGT.m'
disp '   - precomputationForFDGT.m'
disp ' '
disp 'It should be easier to open and run the demo files first because they illustrate'
disp 'how to use the functions. For explanation of each function, please open the file'
disp 'or type "help functionName" into the command line (for example "help DGT").'
disp ' '
disp ' '
disp '3. EXECUTION IN CODE OCEAN'
disp ' '
disp 'In Code Ocean, all figures are saved by the special function codeOceanFigSave.m.'
disp 'In your own MATLAB, these figures are not saved but appear in figure windows as usual.'
disp ' '
disp ' '
disp ' '
disp '--------------------------------------------------------------------------------'
disp ' '
disp ' '


%% run "demo1_DGTusage.m"
disp 'Running "demo1_DGTusage.m" ...'
disp ' '
disp '   - Spectrogram of a speech signal is calculated.'
disp '   - Figure 1 (Spectrogram of the signal) is plotted.'
disp '   - Inverse transform from spectrogram to signal is calculated.'
disp '   - Reconstruction error of the inverse transform is calculated.'
disp '   - Figure 2 (Original and reconstructed signal with their difference) is plotted.'
disp ' '

cd DGT+iPC
run demo1_DGTusage

disp ' '
disp ' '
disp '--------------------------------------------------------------------------------'
disp ' '
disp ' '


%% run "demo2_windowUsage.m"
disp 'Running "demo2_windowUsage.m" ...'
disp ' '
disp '   - A cosine-series-based window is generated.'
disp '   - The corresponding canonical dual window is calculated.'
disp '   - The corresponding canonical tight window is calculated.'
disp '   - Figure 1 (Canonical dual and tight window of the original window) is plotted.'
disp '   - Reconstruction errors for four combination of the windows are calculated.'
disp ' '

run demo2_windowUsage

disp ' '
disp ' '
disp '--------------------------------------------------------------------------------'
disp ' '
disp ' '


%% run "demo3_IFandIPC.m"
disp 'Running "demo3_IFandIPC.m" ...'
disp ' '
disp '   - A cosine-series-based window and its derivative is generated.'
disp '   - Figure 1 (Generated window and its derivative) is plotted.'
disp '   - A numerical derivative of the window is calculated.'
disp '   - Figure 2 (Comparison of analytic and numerical derivatives) is plotted.'
disp '   - Spectrograms are calculated by the window and its derivative.'
disp '   - Bin-wise instantaneous frequency is calculated using the spectrograms.'
disp '   - Instantaneous-phase-corrected (iPC) spectrogram [1] is calculated.'
disp '   - Figure 3 (Calculated spectrograms) is plotted.'
disp '   - Inverse transform from iPC spectrogram to signal is calculated.'
disp '   - Reconstruction error of the inverse transform is calculated.'
disp ' '

run demo3_IFandIPC

disp ' '
disp ' '
disp '--------------------------------------------------------------------------------'
disp ' '
disp ' '


%% run "demo_FDGTusage.m"
disp 'Running "demo_FDGTusage.m" ...'
disp ' '
disp '   - Precomputation necessary for FDGT and invFDGT is performed.'
disp '   - FDGT and invFDGT are repeatedly applied inside loop.'
disp '   - DGT and invDGT are repeatedly applied inside loop for comparison.'
disp ' '

cd fasterDGT
run demo_FDGTusage

disp ' '
disp ' '
disp '--------------------------------------------------------------------------------'
disp ' '
disp ' '

cd ../../