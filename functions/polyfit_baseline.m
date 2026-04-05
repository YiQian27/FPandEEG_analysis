function [detrend_465, detrend_465_2] = polyfit_baseline (signal_405_1,signal_465_1,signal_405_2,signal_465_2,signal_fs)
MeanFilterOrder = 1000; % for smoothing
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;
mouse{14}=1:length(signal_465_1)./signal_fs;

reg = polyfit(signal_405_1(round(mouse{14}*signal_fs)), signal_465_1(round(mouse{14}*signal_fs)), 1);
a = reg(1);
b = reg(2);
    
controlFit = a.*signal_405_1 + b;
controlFit =  filtfilt(MeanFilter,1,double(controlFit));
normDat = (signal_465_1 - controlFit)./controlFit;
delta_465 = normDat * 100;

fs_signal = 1:1:length(signal_465_1);
sec_signal = fs_signal/signal_fs; % time vector for fiber photometry signal

figure
a = subplot(4,1,1);
plot(sec_signal(1000:end), signal_405_1(1000:end));
title('raw control');
b = subplot(4,1,2);
plot(sec_signal(1000:end), signal_465_1(1000:end));
title('raw signal');
c = subplot(4,1,3);
plot(sec_signal(1000:end), signal_465_1(1000:end));
hold on
plot(sec_signal(1000:end), controlFit(1000:end));
title('fitted control');
d = subplot(4,1,4);
plot(sec_signal(1000:end), delta_465(1000:end));
title('normalized signal');
linkaxes([a,b,c,d],'x'); 

ds_factor_FP = 100; % also used for plotting later (section 9b)

% downsampling traces for plotting
delta465_filt = filtfilt(MeanFilter,1,double(delta_465));
ds_delta465_filt = downsample(delta465_filt, ds_factor_FP);
ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting
detrend_465 =detrend(ds_delta465_filt,2,12*3600*fs_signal,"SamplePoints",ds_sec_signal,"Continuous",false);
figure;
plot(detrend_465);

mouse{14}=1:length(signal_465_2)/signal_fs;
reg = polyfit(signal_405_2(round(mouse{14}*signal_fs)), signal_465_2(round(mouse{14}*signal_fs)), 1);
a = reg(1);
b = reg(2);
controlFit_2 = a.*signal_405_2 + b;
controlFit_2 =  filtfilt(MeanFilter,1,double(controlFit_2));
normDat_2 = (signal_465_2 - controlFit_2)./controlFit_2;
delta_465_2 = normDat_2 * 100;

figure
a = subplot(4,1,1);
plot(sec_signal(1000:end), signal_405_2(1000:end));
title('raw control');
b = subplot(4,1,2);
plot(sec_signal(1000:end), signal_465_2(1000:end));
title('raw signal');
c = subplot(4,1,3);
plot(sec_signal(1000:end), signal_465_2(1000:end));
hold on
plot(sec_signal(1000:end), controlFit_2(1000:end));
title('fitted control');
d = subplot(4,1,4);
plot(sec_signal(1000:end), delta_465_2(1000:end));
title('normalized signal');
linkaxes([a,b,c,d],'x'); 

% downsampling traces for plotting
delta465_filt_2 = filtfilt(MeanFilter,1,double(delta_465_2));
ds_delta465_filt_2 = downsample(delta465_filt_2, ds_factor_FP);
ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting
detrend_465_2 =detrend(ds_delta465_filt_2,2,12*3600*fs_signal,"SamplePoints",ds_sec_signal,"Continuous",false);
figure;
plot(detrend_465_2);
