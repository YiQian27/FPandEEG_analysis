function plot_sleep_data(ds_delta465_filt_1, ds_delta465_filt_2, EMG, EEG, wake_woMA_binary_vector, sws_binary_vector, REM_binary_vector, MA_binary_vector, EEG_fs,EMG_fs,ds_signal_fs)

sleepscore_time = 0:length(wake_woMA_binary_vector)-1; % should be same length for wake/sws/REM
fs_signal = 1:1:length(ds_delta465_filt_1);
ds_sec_signal = fs_signal/ds_signal_fs; % time vector for fiber photometry signal
EEG_time = (0:length(EEG)-1)/EEG_fs;
EMG_time = (0:length(EMG)-1)/EMG_fs;
fig = figure;
a = subplot(4,1,1);
    plot_sleep(ds_sec_signal, ds_delta465_filt_1 , sleepscore_time, wake_woMA_binary_vector, sws_binary_vector, REM_binary_vector, MA_binary_vector);
    title('dF/F');
    ylim([-5 5])
b = subplot(4,1,2);
    plot_sleep(ds_sec_signal, ds_delta465_filt_2, sleepscore_time, wake_woMA_binary_vector, sws_binary_vector, REM_binary_vector, MA_binary_vector);
    title('dF/F');
    ylim([-5 5])
c = subplot(4,1,3);
    ds_EMG_time = downsample(EMG_time, 10);
    ds_EMG_rawtrace = downsample(EMG, 10);
    plot_sleep(ds_EMG_time, ds_EMG_rawtrace, sleepscore_time, wake_woMA_binary_vector, sws_binary_vector, REM_binary_vector, MA_binary_vector);
    xlabel('time (s)');
    ylabel('EMG (V)');
d = subplot(4,1,4);
    ds_EEG_rawtrace = downsample(EEG, 10);
    ds_EEG_time = downsample(EEG_time, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time, wake_woMA_binary_vector, sws_binary_vector, REM_binary_vector, MA_binary_vector);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c,d],'x');
