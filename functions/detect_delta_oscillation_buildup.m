function [ result_output ] = detect_delta_oscillation_buildup( ds_delta465_filt_1, fs, NREMinclMA_periods, wake_woMA_periods, defined_interval)
%% Band-pass filter (0.01–0.02 Hz)
signal = ds_delta465_filt_1;
lowcut = 0.01;
highcut = 0.02;
order = 3;
nyquist = fs / 2;
[b, a] = butter(order, [lowcut highcut] / nyquist, 'bandpass');
filtered_signal = filtfilt(b, a, signal);
analytic_signal = hilbert(filtered_signal);
envelope = abs(analytic_signal);
envelope_mean = movmean(envelope, 50 * fs);

%% Collect envelope during NREMinclMA
NREMinclMA_onset = NREMinclMA_periods(:,1);
NREMinclMA_offset = NREMinclMA_periods(:,2);
NE_NREMinclMA = [];

for i = 1:size(NREMinclMA_periods, 1)
    t1 = floor((NREMinclMA_onset + 1) * fs);
    t2 = floor(NREMinclMA_offset * fs);

    if t2 > length(envelope_mean)
        t2 = length(envelope_mean);
    end

    NE_NREMinclMA = [NE_NREMinclMA, envelope_mean(t1:t2)];
end

%% Thresholds
threshold_envelop_lower = prctile(NE_NREMinclMA, 10);
threshold_envelop_upper = prctile(NE_NREMinclMA, 30);

[indices_below_threshold_onset, indices_below_threshold_offset] = ...
    binary_to_OnOff(envelope_mean < threshold_envelop_lower);

[indices_upper_threshold_onset, indices_upper_threshold_offset] = ...
    binary_to_OnOff(envelope_mean > threshold_envelop_upper);

below_threshold_index = find(envelope_mean < threshold_envelop_lower);

%% Filter onsets/offsets within NREM
filtered_offsets = [];
filtered_onsets = [];

for n = 1:length(NREMinclMA_onset)

    range_start = NREMinclMA_onset(n) * fs;
    range_end   = NREMinclMA_offset(n) * fs;

    offsets_in_range = indices_below_threshold_offset( ...
        indices_below_threshold_offset >= range_start & ...
        indices_below_threshold_offset <= range_end);

    first_low_in_NREM = below_threshold_index( ...
        below_threshold_index >= range_start & ...
        below_threshold_index <= range_end);

    if ~isempty(offsets_in_range)
        filtered_offsets = [filtered_offsets, offsets_in_range(1)];
    end

    if ~isempty(first_low_in_NREM)
        filtered_onsets = [filtered_onsets, first_low_in_NREM(1)];
    end
end

%% oscillation detection
Ture_oscllation_start_index = [];
Ture_oscllation_end_index   = [];
Period_osci_builup          = [];
fail_buildup                = [];
Period_wake_before          = [];
Period_osci_before_osci     = [];
Period_osci_after_osci      = [];
NREMinclMA_no_minimum       = [];

for i = 1:size(wake_woMA_periods, 1)

    wake_end = wake_woMA_periods(i, 2);
    NREM_idx = find(NREMinclMA_onset == wake_end);

    if isempty(NREM_idx)
        continue
    end

    start_NREM_index = NREMinclMA_onset(NREM_idx);
    end_NREM_index   = NREMinclMA_offset(NREM_idx);
    % the lowest point should be in first min
    start_index = wake_end * fs;
    end_index   = (wake_end + 60) * fs;

    relevant_onsets = filtered_onsets( ...
        filtered_onsets >= start_index & ...
        filtered_onsets <= end_index);

    if isempty(relevant_onsets)
        NREMinclMA_no_minimum = [NREMinclMA_no_minimum;
                                 start_NREM_index, end_NREM_index];
        continue
    end

    offset_idx = relevant_onsets(1);

    oscillation_starts = find(indices_upper_threshold_onset > offset_idx);
    found_valid_oscillation = false;

    for k = 1:length(oscillation_starts)

        osc_start = indices_upper_threshold_onset(oscillation_starts(k));
        osc_end   = indices_upper_threshold_offset(oscillation_starts(k));
        duration  = osc_end - osc_start;

        if duration > 50 * fs && (osc_start / fs) < end_NREM_index

            Ture_oscllation_start_index = [Ture_oscllation_start_index, osc_start];
            Ture_oscllation_end_index = [Ture_oscllation_end_index, osc_end];
            % the unconsolidated dynamics in mixture
            Period_osci_builup = [Period_osci_builup, (osc_start - start_index) / fs];
            Period_osci_before_osci =  [Period_osci_before_osci; start_index, osc_start];
            % the consolidated dynamics in mixture
            Period_osci_after_osci = [Period_osci_after_osci; osc_start, end_NREM_index * fs];

            Period_wake_before = [Period_wake_before; wake_woMA_periods(i,:)];
            found_valid_oscillation = true;
            break
        end
    end

    if ~found_valid_oscillation
        fail_buildup = [fail_buildup; start_NREM_index, end_NREM_index];
    end
end


time_in_seconds = Ture_oscllation_start_index ./ fs;
interval = defined_interval * 3600;  % Change!

first_indices_buildup = Period_osci_builup(time_in_seconds <= interval);
second_indices_buildup = Period_osci_builup(time_in_seconds > interval);
if ~isempty(fail_buildup)
    time_in_seconds_2 = fail_buildup(:,1);
    first_indices_fail = fail_buildup(time_in_seconds_2 <= interval,:);
    second_indices_fail = fail_buildup(time_in_seconds_2 > interval,:);

        
    Times_fail_1 = size(first_indices_fail,1);
    Times_fail_2 = size(second_indices_fail,1);
        if ~isempty(first_indices_fail)
            Mean_NREMduration_failure_1 = mean(first_indices_fail(:,2)-first_indices_fail(:,1));
        else 
            Mean_NREMduration_failure_1 = NaN;
        end

        if ~isempty(second_indices_fail)
            Mean_NREMduration_failure_2 = mean(second_indices_fail(:,2)-second_indices_fail(:,1));
        else 
            Mean_NREMduration_failure_2 = NaN;
        end
else
    Times_fail_1 = 0;
    Times_fail_2 = 0;
    Mean_NREMduration_failure_2 = NaN;
    Mean_NREMduration_failure_1 = NaN;
end




time_in_seconds_no_min = NREMinclMA_no_minimum(:,1);
first_no_min = NREMinclMA_no_minimum(  time_in_seconds_no_min <= interval, :);

second_no_min = NREMinclMA_no_minimum( time_in_seconds_no_min > interval, :);



Times_buildup_1 = length(first_indices_buildup);
Times_buildup_2 = length(second_indices_buildup);

Mean_oscillation_buildup_1 = mean (first_indices_buildup);
Mean_oscillation_buildup_2 = mean (second_indices_buildup);
result_output=[Times_fail_1, Mean_NREMduration_failure_1, Times_buildup_1, Mean_oscillation_buildup_1,size(first_no_min,1),Times_fail_2, Mean_NREMduration_failure_2, Times_buildup_2, Mean_oscillation_buildup_2,size(second_no_min,1)];

end