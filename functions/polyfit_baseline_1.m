function [ds_sec_signal,ds_detrend_465] = polyfit_baseline_1 (signal_405_1,signal_465_1,signal_fs)
MeanFilterOrder = 1000; % for smoothing
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;
mouse{14}=1:length(signal_465_1)/signal_fs;

reg = polyfit(signal_405_1(round(mouse{14}*signal_fs)), signal_465_1(round(mouse{14}*signal_fs)), 1);
a = reg(1);
b = reg(2);
    
controlFit = a.*signal_405_1 + b;
controlFit =  filtfilt(MeanFilter,1,double(controlFit));
normDat = (signal_465_1 - controlFit)./controlFit;
delta_465 = normDat * 100;

fs_signal = 1:1:length(signal_465_1);
sec_signal = fs_signal/signal_fs; % time vector for fiber photometry signal
% 
% figure
% a = subplot(4,1,1);
% plot(sec_signal(1000:end), signal_405_1(1000:end));
% title('raw control');
% b = subplot(4,1,2);
% plot(sec_signal(1000:end), signal_465_1(1000:end));
% title('raw signal');
% c = subplot(4,1,3);
% plot(sec_signal(1000:end), signal_465_1(1000:end));
% hold on
% plot(sec_signal(1000:end), controlFit(1000:end));
% title('fitted control');
% d = subplot(4,1,4);
% plot(sec_signal(1000:end), delta_465(1000:end));
% title('normalized signal');
% linkaxes([a,b,c,d],'x'); 

% downsampling traces for plotting
ds_factor_FP =100;
delta465_filt = filtfilt(MeanFilter,1,double(delta_465));
ds_delta465 = downsample(delta465_filt, ds_factor_FP);
ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting
ds_detrend_465 =detrend(ds_delta465,2,12*3600*fs_signal,"SamplePoints",ds_sec_signal,"Continuous",false);
% figure;
% plot(ds_detrend_465);
