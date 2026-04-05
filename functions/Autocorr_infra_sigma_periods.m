function [final_period,period_distribution] = Autocorr_infra_sigma_periods(signal_trace, fs_ref, analysis_periods, varargin)
% estimate_signal_period estimates the main oscillatory period of a signal based on autocorrelation.
%
% INPUTS:
%   signal_trace         - 1D array of signal data
%   fs_ref               - frequency fot the signal_trace
%   analysis_periods     - [n x 2] matrix, start and end times (in seconds) for periods to analyze
%
% OPTIONAL ARGUMENTS:
%   'MinPeriodDur'       - Minimum valid segment duration (seconds) [default = 120]
%   'MaxLagTime'         - Maximum lag time considered in autocorrelation (seconds) [default = 120]
%   'HarmonicsTolerance' - Tolerance for harmonic rejection [default = 0.3]
%   'MinPeakProminence'  - Minimum prominence of peaks [default = 0.15]
%   'SmoothingFactor'    - Smoothing factor for autocorrelation [default = 1.5]
%   'MinPeakDistance'    - Minimum distance between peaks (seconds) [default = 1]
%   'MinPeakHeight'      - Minimum peak height [default = 0]
%   'PlotIndividual'     - Plot each segment's autocorrelation and period [default = true]
%
% OUTPUT:
%   final_period         - Weighted average estimated period across all segments (in seconds)

%% Parse inputs
p = inputParser;
addParameter(p, 'MinPeriodDur', 120);
addParameter(p, 'MaxLagTime', 120);
addParameter(p, 'HarmonicsTolerance', 0.3); % test and change
addParameter(p, 'MinPeakProminence', 0.05);
addParameter(p, 'SmoothingFactor', 5);
addParameter(p, 'MinPeakDistance', 1);
addParameter(p, 'MinPeakHeight', 0);
addParameter(p, 'PlotIndividual', true);
parse(p, varargin{:});
params = p.Results;

%% Setup basic parameters
t1 = analysis_periods(:,1); % Start times
t2 = analysis_periods(:,2); % End times
tsamp1 = floor(t1 * fs_ref); % Start indices
tsamp2 = floor(t2 * fs_ref); % End indices

max_lag_points = floor(params.MaxLagTime * fs_ref); % Max lags for xcorr
lags_common = -max_lag_points:max_lag_points; % Common lags across segments

% Initialize storage
autocorr_results = {}; 
period_durations = [];
harmonic_ratios = [];
master_periods = [];
period_weights = [];

%% Loop through each analysis period
for i = 1:numel(tsamp1)
    if (tsamp2(i) - tsamp1(i)) < params.MinPeriodDur * fs_ref
        continue;
    end

    tsamp1(i) = max(tsamp1(i), 1);
    tsamp2(i) = min(tsamp2(i), length(signal_trace));
    % detrend
    current_signal = signal_trace(tsamp1(i):tsamp2(i));
    [pfit, ~, mu] = polyfit((1:numel(current_signal))', current_signal, 5);  % try 3?
    detrended = current_signal - polyval(pfit, (1:numel(current_signal))', [], mu)';
    
    % Calculate autocorrelation
    [autocorr_result, ~] = xcorr(detrended, max_lag_points, 'coeff');
    % Store results
    autocorr_results{end+1} = autocorr_result;
    period_durations(end+1) = (tsamp2(i) - tsamp1(i)) / fs_ref;
end

%% Analyze autocorrelation results
for segIdx = 1:length(autocorr_results)
    if isempty(autocorr_results{segIdx}) || period_durations(segIdx) < params.MinPeriodDur
        continue;
    end
    
    raw_autocorr = autocorr_results{segIdx};
    smoothed_autocorr = imgaussfilt(raw_autocorr, params.SmoothingFactor);
    time_lags = lags_common / fs_ref;
    
    % Find peaks
    [peaks, locs] = findpeaks(smoothed_autocorr, time_lags, ...
        'MinPeakProminence', params.MinPeakProminence, ...
        'MinPeakDistance', params.MinPeakDistance, ...
        'MinPeakHeight', params.MinPeakHeight);
    
    % symmetry, positive half
    valid = locs > 10 & peaks > 0;
    peaks = peaks(valid);
    locs = locs(valid);
    
    if isempty(locs)
        continue;
    end
    
    % Use the smallest positive lag as base period
    [base_period, base_idx] = min(locs);
    candidate_periods = locs;
    candidate_weights = peaks / max(peaks); % for weighting peaks within one peiod 
    
    % Identify and remove harmonics
    is_harmonic = false(size(candidate_periods));
    for j = 1:length(candidate_periods)
        if j == base_idx
            continue;
        end
        ratio = candidate_periods(j) / base_period;
        nearest_int = round(ratio); % closest periods
        if abs(candidate_periods(j)/nearest_int - base_period)/base_period <= params.HarmonicsTolerance
            is_harmonic(j) = true;
            harmonic_ratios = [harmonic_ratios; nearest_int];
        end
    end
    
    valid_periods = candidate_periods(~is_harmonic);
    valid_weights = candidate_weights(~is_harmonic);
    
    % Calculate weighted period for this segment
    if ~isempty(valid_periods)
        weighted_period = sum(valid_periods .* valid_weights) / sum(valid_weights); % first weighting -- within one NREMinclMA period
        master_periods(end+1) = weighted_period;
        period_weights(end+1) = period_durations(segIdx);
    end
    
    % Plot individual results
    if params.PlotIndividual
        figure('Position', [200 200 800 400], 'Color','w');
        
        % Plot autocorrelation
        subplot(1,2,1);
        plot(time_lags, smoothed_autocorr, 'LineWidth', 1.5, 'Color', [0.2 0.6 0.8]);
        hold on;
        scatter(locs, peaks, 40, 'r', 'filled');
        title(sprintf('Segment %d Autocorrelation', segIdx));
        xlabel('Lag Time (s)');
        ylabel('Autocorrelation Coefficient');
        grid on;
        
        % Plot periods
        subplot(1,2,2);
        if ~isempty(valid_periods)
            bar(valid_periods, valid_weights, 'FaceColor', [0.8 0.4 0.1]);
            hold on;
            plot(weighted_period*[1 1], ylim(), '--r', 'LineWidth', 2);
            title(sprintf('Weighted Period: %.2f s', weighted_period));
            xlabel('Period (s)');
            ylabel('Weight');
        else
            text(0.5, 0.5, 'No valid period', 'HorizontalAlignment', 'center');
            title('Period Results');
        end
        set(gca, 'XLim', [0 max(10, 2*base_period)]);
        grid on;
    end
end

%% Final calculation
if ~isempty(master_periods)
    final_period = sum(master_periods .* period_weights) / sum(period_weights); % second weighting -- between NREMinclMA periods
    fprintf('Final weighted average period: %.2f seconds\n', final_period);
    fprintf('Harmonic ratios distribution:\n');
    tabulate(harmonic_ratios);
else
    error('No effective periods detected.');
end

figure;
histogram(master_periods, 'BinWidth', 10);

h = histogram(master_periods, 'BinWidth', 10);

bin_edges = h.BinEdges;          
bin_starts = bin_edges(1:end-1); 
bin_counts = h.Values;
period_distribution = [bin_starts(:), bin_counts(:)];


end