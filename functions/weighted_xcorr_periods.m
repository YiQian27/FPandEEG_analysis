%% cross_correlation 
function [mean_cc_weighted, lags_sec, lag_at_max, cc1_clean, lags_clean] = weighted_xcorr_periods(signal1, frequency1, signal2, frequency2, period, TimeLagInSamples)
% TimeLagInSamples = 1000; 
% period = NREMinclMA_periods_day1;
% signal1=ds_delta465_filt_2_day1;
% signal2=ds_delta465_filt_1_day1;
% frequency1 = fs_ref_NE ; 
% frequency2 = fs_ref_NE ; 
cc1 = [];
lags = [];
Cor_NREMinclMA_periods =[];
for i = 1:length(period)
    start_time_1 = round(period(i,1)*frequency1+1);
    stop_time_1 = round(period(i,2)*frequency1);
    if stop_time_1 > length(signal1)
        stop_time_1 = length(signal1);
    end
    period_NREM_signal1 = signal1(1, start_time_1:stop_time_1);

    start_time_2 = round(period(i,1)*frequency2+1);
    stop_time_2 = round(period(i,2)*frequency2);
    if stop_time_2 > length(signal2)
        stop_time_2 = length(signal2);
    end
    period_NREM_signal2 = signal2(1, start_time_2:stop_time_2);

    target_len = length(period_NREM_signal1);
    orig_len = length(period_NREM_signal2);

    x = 1:orig_len;
    xq = linspace(1, orig_len, target_len);

    reampled_NE = interp1(x, period_NREM_signal2, xq, 'linear');

    [cc, lags_vec] = xcorr(unity(period_NREM_signal1), unity(reampled_NE), round(TimeLagInSamples), 'unbiased');
    Cor_NREMinclMA_periods = [Cor_NREMinclMA_periods, period(i,2)-period(i,1)];
    cc1(i,:) = cc;
    lags(i,:) = lags_vec;

end

nan_rows = any(isnan(cc1), 2);  

cc1_clean = cc1(~nan_rows, :);   
lags_clean = lags(~nan_rows,:);
Cor_NREMinclMA_periods_clean = Cor_NREMinclMA_periods(~nan_rows); 

weights = Cor_NREMinclMA_periods_clean(:);             
weights = weights / sum(weights);                       

mean_cc_weighted = weights' * cc1_clean;               

%plot
figure;
plot(lags(1,:)./frequency1, mean_cc_weighted, 'k', 'LineWidth', 2);
xlabel('Lag (s)');
ylabel('Weighted Mean Cross-correlation');
title('Duration-Weighted Mean xcorr (NaN removed)');
lags_sec = lags_clean(1,:) ./ frequency1;
[maxValue, idx] = max(mean_cc_weighted);
lag_at_max = lags(1,idx)./frequency1; 