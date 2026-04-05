%% define mouse data
clear all
close all
% data structure:
    % 1) Mouse data address
    % 2) Mouse sleepscore data address
    % 3) Mouse FP 405 name
    % 4) Mouse FP 465 name
    % 6) recording start time
    % 9) Mouse EEG name
    % 10) Mouse EEG channel
    % 10) Mouse EMG name
Example_ms = {'/Users/qy/Desktop/Master project/NE oscillation/Data/Young/M10_1_783_2_24hr' '/Users/qy/Desktop/Master project/3. NE oscillation/Data/Young/M10_1_783_2_24hr/M10_24hr.exp' 'x65A' 'x05A' 'x' '01-Jan-2024 6:58:39' 0 'x' 'EEGw' 1 'EMG1' 'x' 'x'};

mouse= Example_ms; % Change to the subject you want to analyze
t1={(mouse{3})}; % Time when recording is started
t2={'01-Jan-2024 19:00:00'}; % Change to the day and time you want to start analyzing from
analysis_hours =24; % analysis period
%% load data
data_FPrig = TDTbin2mat(mouse{1});
%% extract channels
signal_fs = data_FPrig.streams.(mouse{3}).fs; % sampling frequency for fiber photometry signal
signal_465= data_FPrig.streams.(mouse{3}).data; %signal
signal_405= data_FPrig.streams.(mouse{4}).data; %isosbetstic control
EEG_fs =  data_FPrig.streams.(mouse{9}).fs; %sampling frequency for EEG signal 

EEG_rawtrace = data_FPrig.streams.(mouse{9}).data; %EEG signal
EEG_rawtrace = EEG_rawtrace(mouse{10},:); %add channel 
EMG_rawtrace =  data_FPrig.streams.(mouse{11}).data; %EMG 
%% time sequences

fs_signal = 1:1:length(signal_465);
sec_signal = fs_signal /signal_fs; % time vector for fiber photometry signal

fs_signal_EEG = 1:1:length(EEG_rawtrace);
sec_signal_EEG = fs_signal_EEG/EEG_fs; % time vector for EEG signal

fs_signal_EMG = 1:1:length(EMG_rawtrace);
sec_signal_EMG = fs_signal_EMG/EEG_fs; % time vector for EMG signal
%% preliminary processing
[ds_sec_signal,ds_detrend_465] = polyfit_baseline_1 (signal_405,signal_465,signal_fs);
%% Time alignment with the priod you interest
t1={(mouse{6})}; % Time when recording is started
t2={'02-Jan-2024 7:04:00'}; % Time for analyzing start
analysis_hours =12; % Analysis hour

interval = datetime(t2)-datetime(t1); 
interval_s = seconds(interval);
time_correction = 0; % change it if you need to correct time for alignment
starttime=interval_s+1+time_correction;
endtime=analysis_hours*3600+interval_s+time_correction;
ds_signal_fs=signal_fs/100;
EEG=EEG_rawtrace(((starttime-1)*EEG_fs+1):endtime*EEG_fs);
EMG=EMG_rawtrace(((starttime-1)*EEG_fs+1):endtime*EEG_fs);
ds_delta465_filt_1=ds_detrend_465(starttime*ds_signal_fs:endtime*ds_signal_fs);

ds_sec_signal_1=1:1:length(ds_delta465_filt_1);
ds_sec_signal_1=ds_sec_signal_1./ds_signal_fs;
fs_signal_EEG_cut = 1:1:length(EEG);
sec_signal_EEG_cut = fs_signal_EEG_cut/EEG_fs; % time vector for EEG signal
frq = EEG_fs; 

%% sleepscore data

ViewpointData.FileInfo=loadEXP([mouse{2}],'no'); %   
TimeReldebSec=0; 

TimeRelEndSec=ViewpointData.FileInfo.BinFiles.Duration; 

[Data,Time]=ExtractContinuousData([],ViewpointData.FileInfo,[],TimeReldebSec, TimeRelEndSec,[],1);
[FullHypno,TimeScaleAbs,TimeScaleBin,TimeScaleHypno]=ExtractFullHypno(ViewpointData,1);

starttime = interval_s+1+time_correction;
endtime = analysis_hours*3600+interval_s+time_correction;
FullHypno = FullHypno(starttime:endtime);
%% Sleep composition
[Sleep_Composition, wake_binary_vector, sws_binary_vector, REM_binary_vector, NREMinclMA_binary_vector,MA_binary_vector,wake_woMA_binary_vector] = sleepcompostion_computing(FullHypno);

[wake_onset, wake_offset] = binary_to_OnOff(wake_binary_vector);
[sws_onset, sws_offset] = binary_to_OnOff(sws_binary_vector);
[REM_onset, REM_offset] = binary_to_OnOff(REM_binary_vector);
[wake_woMA_onset, wake_woMA_offset] = binary_to_OnOff(wake_woMA_binary_vector);
[NREMinclMA_onset, NREMinclMA_offset] = binary_to_OnOff(NREMinclMA_binary_vector);


wake_periods = [wake_onset wake_onset+wake_duration];
sws_periods = [sws_onset sws_onset+sws_duration];
REM_periods = [REM_onset REM_onset+REM_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_offset];
NREMinclMA_periods = [NREMinclMA_onset NREMinclMA_offset];


%% State transitions
ds_delta465_rebuilt = Decom_wav(ds_delta465_filt_1, 6, false);

Period_NREM_for_thresholding = [];
for i = 1:length(NREMinclMA_periods)
    start_time_NE = round(NREMinclMA_periods(i,1) * ds_signal_fs + 1);
    stop_time_NE = round(NREMinclMA_periods(i,2) * ds_signal_fs);
    stop_time_NE = min(stop_time_NE, length(ds_delta465_rebuilt));
    Period_NREM_for_thresholding = [Period_NREM_for_thresholding, ...
                                    ds_delta465_rebuilt(start_time_NE:stop_time_NE)];
end
% min
range_NE = prctile(Period_NREM_for_thresholding(100:end), 90) - prctile(Period_NREM_for_thresholding(100:end), 10);
MinSeparation_min = 15 * ds_signal_fs;   
[inverted_peaks, min_locs] = findpeaks(-ds_delta465_rebuilt, 'MinPeakProminence',0.1*range_NE, 'MinPeakDistance', MinSeparation_min);
min_values = -inverted_peaks;
% max
MinSeparation_peak = 15 * ds_signal_fs;  
[peaks, peak_locs] = findpeaks(ds_delta465_rebuilt, 'MinPeakProminence',0.1*range_NE, 'MinPeakDistance', MinSeparation_peak);

% % For example: NE trace from NREM to MA 
% FullHypno_MA=FullHypno;
% for i = 1:length(MA_periods)
%     FullHypno_MA((MA_periods(i,1)+1):MA_periods(i,2))=3;
% end
% 
% turning_sws_MA=[];
% for i = 1:length(FullHypno_MA)-1
%     if  FullHypno_MA(i) == 2 &&  FullHypno_MA(i+1) == 3
%         turning_sws_MA = [turning_sws_MA, i + 1]; 
%     end
% end
% 
% deNE_fre = ds_signal_fs;
% MA_onset_in_NE = round(turning_sws_MA.* deNE_fre);
% last_min_before_MA = zeros(length(MA_onset_in_NE), 1);
% for i = 1:length(MA_onset_in_NE)
%     min_before_MA = min_locs(min_locs < MA_onset_in_NE(i));
%     if ~isempty(min_before_MA) & min_before_MA(end)- MA_onset_in_NE(i) < 30 * deNE_fre
%         last_min_before_MA(i) = min_before_MA(end);
%     else
%         last_min_before_MA(i) = NaN;
%     end
% end
% 
% fs = ds_signal_fs;          
% search_win_sec = 0.5;         
% search_win = round(search_win_sec * fs);
% raw_signal = ds_delta465_filt_1;
% N = length(raw_signal);
% 
% last_min_before_MA_raw_vals = nan(size(last_min_before_MA));
% last_min_before_MA_raw_locs = nan(size(last_min_before_MA));
% 
% for i = 1:length(last_min_before_MA)
%     if isnan(last_min_before_MA(i))
%         continue
%     end
% 
%     center = last_min_before_MA(i);
%     idx = max(1, center-search_win) : min(N, center+search_win);
%     [last_min_before_MA_raw_vals(i), idx_rel] = min(raw_signal(idx));
%     last_min_before_MA_raw_locs(i) = idx(idx_rel);
% end
% 
% st = round(50 * deNE_fre);
% sp = round(75 * deNE_fre);
% 
% NE_sws_MA=[];
% for i=1:length(last_min_before_MA_raw_locs)
%     if last_min_before_MA_raw_locs(i)-st>0 && last_min_before_MA_raw_locs(i)+sp<length(ds_delta465_filt_1)
%         NE_sws_MA=[NE_sws_MA;ds_delta465_filt_1(last_min_before_MA_raw_locs(i)-st:last_min_before_MA_raw_locs(i)+sp)];
%     end
% end
% 
% NE_sws_MA_mean=mean(NE_sws_MA,1);
% NE_sws_MA_mean=NE_sws_MA_mean';
% amplitude_NREM_MA = max(NE_sws_MA_mean) - min(NE_sws_MA_mean);
% Nor_NE_sws_MA_mean = NE_sws_MA_mean - mean(NE_sws_MA_mean(1:round(25 * ds_signal_fs)));
% 
% sec_NE_event=[1:1:length(Nor_NE_sws_MA_mean)]./deNE_fre;
% figure;
% plot(sec_NE_event,Nor_NE_sws_MA_mean,'b');
% hold on;
% plot([50 50],[-3 3],'--')
% hold off;

%% Increasing Slope 
[t_slope, slope] = derivative (ds_delta465_filt_1, ds_signal_fs, 4, 1, true);

local_max_indices = islocalmax(slope);
local_max_slope = slope(local_max_indices);

Positive_local_max_indices = local_max_indices & (slope > prctile(slope,80));
Positive_local_max_slope = slope(Positive_local_max_indices);

%% infraslow NE spectrum and periods
[infra_fre_NE, spectrum_epoch_weighted_NE] = analyze_ISO_spectrum(ds_delta465_filt_1, ds_signal_fs, NREMinclMA_periods, 120, 5, true);
final_period_NE = Autocorr_infra_periods(ds_delta465_filt_1, ds_signal_fs, NREMinclMA_periods);
% Continued infraslow oscillation power
lowcut = 0.01;
highcut = 0.02;
order = 3;
nyquist = ds_signal_fs/2;
[b,a] = butter(order,[lowcut highcut]/nyquist,'bandpass');
filtered_signal = filtfilt(b,a,ds_delta465_filt_1);
