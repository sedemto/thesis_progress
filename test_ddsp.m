addpath(genpath("PEMO-Q_v1.4.1"))
[orig, fsref] = audioread("for_ddsp/input2_original.wav");
nobatch = audioread("for_ddsp/output2_violin_noBatch.wav");
batch16 = audioread("for_ddsp/output2_violin_Batch16.wav");
batch100 = audioread("for_ddsp/output2_violin_Batch100.wav");
batch1_exp = audioread("for_ddsp/output2_violin_noBatch_exp.wav");
batch1_mfccTime =  audioread("for_ddsp/output2_violin_noBatch_mfccTime.wav");


SNR_nobatch = snr(orig, orig-nobatch);
SNR_batch16 = snr(orig, orig-batch16);
SNR_batch100 = snr(orig, orig-batch100);
SNR_nobatch_exp = snr(orig, orig-batch1_exp);
SNR_nobatch_mfccTime = snr(orig, orig-batch1_mfccTime);


orig_44k = resample(orig,44100,fsref);

nobatch_44 = resample(nobatch,44100,fsref);
batch16_44 = resample(batch16,44100,fsref);
batch100_44 = resample(batch100,44100,fsref);
batch1_exp_44 = resample(batch1_exp,44100,fsref);
batch1_mfccTime_44  = resample(batch1_mfccTime,44100,fsref);

[PSM,PSMt,ODG_nobatch,PSM_inst] = audioqual_silent(orig_44k,nobatch_44,44100);
[PSM,PSMt,ODG_batch16,PSM_inst] = audioqual_silent(orig_44k,batch16_44,44100);
[PSM,PSMt,ODG_batch100,PSM_inst] = audioqual_silent(orig_44k,batch100_44,44100);
[PSM,PSMt,ODG_batch1_exp,PSM_inst] = audioqual_silent(orig_44k,batch1_exp,44100);
[PSM,PSMt,ODG_batch1_mfccTime,PSM_inst] = audioqual_silent(orig_44k,batch1_mfccTime,44100);

table = ["input1" "noBatch" "Batch16" "Batch100" "batch1_exp" "batch1_mfccTime";"SNR" SNR_nobatch SNR_batch16 SNR_batch100 SNR_nobatch_exp SNR_nobatch_mfccTime;"ODG" ODG_nobatch ODG_batch16 ODG_batch100 ODG_batch1_exp ODG_batch1_mfccTime];