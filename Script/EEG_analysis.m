%% power
% for example delta
Data_EEG = double(EEG); % EEG trace
delta_band = bandpass(Data_EEG', [1 4], EEG_fs);
delta_power_highfre = abs(hilbert(delta_band)).^2;
resample_fre = 5;
delta_power = resample(delta_power_highfre, resample_fre, round(EEG_fs));
%% PSD
frw = 0:0.1:30;
window = 5; %second
window_samples = round(EEG_fs * window); 
[transition_spectrogram, F, T] = spectrogram(Data_EEG, window_samples, [], frw, frq, 'yaxis');
mean_spectrogram = log(abs(transition_spectrogram));
%% sigma infraslow 
NREMinclMA_periods = NREMinclMA_periods; % corresponding NREMinclMA periods, (NREMinclMA_onset,NREMinclMA_offset), seconds
resample_fre=5;
sigma_band = bandpass(Data_EEG, [10 15], EEG_fs);
signal_trace= abs(hilbert(sigma_band)).^2;
t = (0:length(signal_trace)-1) / EEG_fs;
sec_idx = floor(t*resample_fre) + 1;
sigma_power = accumarray(sec_idx(:), signal_trace(:), [], @mean);

[infra_fre, spectrum_epoch_weighted] = analyze_ISO_spectrum(signal_trace, EEG_fs, NREMinclMA_periods, 120, 5, true);
final_period_sigma = Autocorr_infra_sigma_periods(sigma_power', resample_fre, NREMinclMA_periods);
%% spindle_related
output = detect_spindles_full( ...
    EEG, EEG_fs, ...
    sws_periods;