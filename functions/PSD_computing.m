function [mean_PXX, PXX, f] = PSD_computing(EEG, analysis_period, frq, freq_band)

% analysis_period: [start end] in sleep-score seconds
analysis_period_EEG = [ ...
   floor(analysis_period(:,1) * frq+1), ...
     floor(analysis_period(:,2) * frq) ];
analysis_duration = analysis_period(:,2) - analysis_period(:,1);

PXX = [];

for i = 1:size(analysis_period_EEG,1)

    idx1 = analysis_period_EEG(i,1);
    idx2 = min(analysis_period_EEG(i,2), length(EEG));

    segment = EEG(idx1:idx2);

    [pxx, f] = pwelch(double(segment), [], [], freq_band, frq);
    PXX(:,i) = pxx;
end

weights = analysis_duration ./ sum(analysis_duration);
mean_PXX = sum(PXX .* weights', 2);

end