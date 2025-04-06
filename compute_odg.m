clc;clear all; close all;
addpath(genpath("PEMO-Q_v1.4.1"))


folderToIt = "C:\Users\petob\OneDrive\Documents\SCHOOL&WORK\PHAIN_dp\PHAIN-v1.0.1\TomoroTanaka-PHAIN-e9f0099\results";

NN = 8;
M = 6;
odgResults8 =NaN(NN, M, 1);

for nn=1:NN
    example = "example"+nn;
    reference_path =  "C:\Users\petob\OneDrive\Documents\SCHOOL&WORK\PHAIN_dp\PHAIN-v1.0.1\TomoroTanaka-PHAIN-e9f0099\originals_dpai\audio_original_example"+string(nn-1)+".wav";
    [reference_audio,fsref] = audioread(reference_path);
    ref_44k = resample(reference_audio,44100,fsref);
    for m=1:M
        mask = "mask"+m;
        newFolderPath =  folderToIt+"\"+example+mask+".wav";
        [reconstructed_audio,fsrec] = audioread(newFolderPath);
        rec_44 = resample(reconstructed_audio,44100,fsrec);
        
        [PSM,PSMt,ODG,PSM_inst] = audioqual_silent(ref_44k,rec_44,44100);
        odgResults8(nn,m,1) = ODG;

    end
end
save("resultsPHAIN_ODG_TEST_GCPA_it500_proj_l10e-3.mat","odgResults8");