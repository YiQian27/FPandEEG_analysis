function out = detect_spindles_full( ...
    EEG, EEG_fs, ...
    sws_periods)

% ==============================
% 1. Bandpass filter (spindle band)
% ==============================
fs = EEG_fs;

d = designfilt('bandpassiir', ...
    'StopbandFrequency1', 3/(fs/2), ...
    'PassbandFrequency1', 10/(fs/2), ...
    'PassbandFrequency2', 15/(fs/2), ...
    'StopbandFrequency2', 22/(fs/2), ...
    'StopbandAttenuation1', 24, ...
    'StopbandAttenuation2', 24, ...
    'DesignMethod', 'butter');

filtered_EEG = filtfilt(d, double(EEG));

% ==============================
% 2. RMS envelope
% ==============================
window_size = round(0.75 * fs);
[env, ~] = envelope(filtered_EEG, window_size, 'rms');
cubed_rms = env.^3;

% ==============================
% 3. Extract NREM segments
% ==============================
sws_periods_EEG = round(sws_periods * fs);

filtered_EEG_NREM = arrayfun(@(i) ...
    cubed_rms((sws_periods_EEG(i,1)+1):min(sws_periods_EEG(i,2), length(EEG))), ...
    1:size(sws_periods_EEG,1), 'UniformOutput', false);

NREM_to_EEG_idx = cell(size(sws_periods_EEG,1),1);
for i = 1:size(sws_periods_EEG,1)
    NREM_to_EEG_idx{i} = (sws_periods_EEG(i,1)+1):sws_periods_EEG(i,2);
end

cubed_rms_NREM = cell2mat(filtered_EEG_NREM);
NREM_to_EEG_idx = cell2mat(NREM_to_EEG_idx');

% ==============================
% 4. Sliding window spindle detection
% ==============================
segment_length = round(15 * 60 * fs);
step_size = round(segment_length * 0.5);

segment_starts = 1:step_size:(length(cubed_rms_NREM) - segment_length + 1);

epoch_lengths = cellfun(@length, filtered_EEG_NREM);
epoch_end_idx = cumsum(epoch_lengths);

is_NREM_continuous = true(1, sum(epoch_lengths));
is_NREM_continuous(epoch_end_idx(1:end-1)) = false;

spindle_events_NREM = [];

for s = 1:length(segment_starts)

    seg_start = segment_starts(s);
    seg_end = min(seg_start + segment_length - 1, length(cubed_rms_NREM));

    segment = cubed_rms_NREM(seg_start:seg_end);

    mean_val = mean(segment);
    lower_th = 1.2 * mean_val;
    upper_th = 3.5 * mean_val;

    above_lower = segment > lower_th;
    above_upper = segment > upper_th;

    if ~any(above_upper)
        continue
    end

    [onsets, offsets] = binary_to_OnOff(above_lower);

    for i = 1:length(onsets)

        o = onsets(i)+1;
        f = offsets(i);

        if any(above_upper(o:f))
            event_start = seg_start + o - 1;
            event_end   = seg_start + f - 1;

            if any(~is_NREM_continuous(event_start:event_end-1))
                continue
            end

            spindle_events_NREM = [spindle_events_NREM; event_start, event_end];
        end
    end
end

% ==============================
% 5. Merge + duration filter
% ==============================
spindle_events_NREM = sortrows(spindle_events_NREM);

merged = spindle_events_NREM(1,:);
for i = 2:size(spindle_events_NREM,1)
    curr = spindle_events_NREM(i,:);
    prev = merged(end,:);

    if curr(1) <= prev(2)
        merged(end,2) = max(prev(2), curr(2));
    else
        merged = [merged; curr];
    end
end

spindle_events_NREM = merged;

Duration = (spindle_events_NREM(:,2) - spindle_events_NREM(:,1)) ./ fs;
valid_idx = Duration > 0.5 & Duration < 10;

spindle_events_NREM = spindle_events_NREM(valid_idx,:);

% ==============================
% 6. Map back to EEG index
% ==============================
spindle_events_EEG = zeros(size(spindle_events_NREM));

for i = 1:size(spindle_events_NREM,1)
    spindle_events_EEG(i,1) = NREM_to_EEG_idx(spindle_events_NREM(i,1));
    spindle_events_EEG(i,2) = NREM_to_EEG_idx(spindle_events_NREM(i,2));
end

spindle_on = spindle_events_EEG(:,1);
spindle_off = spindle_events_EEG(:,2);
spindle_center = unique((spindle_on + spindle_off)/2);

% ==============================
% 7. Frequency + amplitude (FFT + RMS)
% ==============================
window_length = 1;
step_size = 0.15;

Nw = round(window_length * fs);
Ns = round(step_size * fs);

signal = filtfilt(d, double(EEG));
num_steps = floor((length(signal) - Nw) / Ns) + 1;

peak_freq = zeros(1, num_steps);
rms_amp = zeros(1, num_steps);

for i = 1:num_steps
    idx = (i-1)*Ns + (1:Nw);
    x = signal(idx);

    rms_amp(i) = rms(x);

    X = abs(fft(x .* hamming(Nw)'));
    f = (0:length(X)-1)*(fs/length(X));

    [~, imax] = max(X);
    peak_freq(i) = f(imax);
end

scale = length(signal) / num_steps;
on_idx = round(spindle_on / scale);
off_idx = round(spindle_off / scale);

num_sp = length(on_idx);

sp_freq = zeros(1,num_sp);
sp_amp  = zeros(1,num_sp);

for i = 1:num_sp
    sp_freq(i) = median(peak_freq(on_idx(i):off_idx(i)));
    sp_amp(i)  = mean(rms_amp(on_idx(i):off_idx(i)));
end

% ==============================
% 8. Local density
% ==============================
local_density = zeros(size(spindle_center));

for i = 1:length(spindle_center)
    win = spindle_center(i) + [-30 30]*fs;
    local_density(i) = sum(spindle_center >= win(1) & spindle_center <= win(2));
end

% ==============================
% 9. Clustering (train detection)
% ==============================
intervals = diff(spindle_on) ./ fs;

clusters = {};
current = [];

for i = 1:length(intervals)
    if intervals(i) <= 6
        if isempty(current)
            current = [i i+1];
        else
            current = [current i+1];
        end
    else
        if ~isempty(current)
            clusters{end+1} = current;
            current = [];
        end
    end
end

cluster_all = cell2mat(clusters);
single_sp = setdiff(1:length(spindle_on), cluster_all);

% ==============================
% 10. Metrics
% ==============================
NumCluster = length(clusters);
AvgPerCluster = mean(cellfun(@numel, clusters));
PropCluster = length(cluster_all)/length(spindle_on)*100;

% ISI
ISI = [];
for i = 1:length(clusters)
    idx = clusters{i};
    ISI = [ISI, intervals(idx(1:end-1))];
end

Mean_ISI = mean(ISI);

% ITI
ITI = [];
for i = 2:length(clusters)
    gap = (spindle_on(clusters{i}(1)) - spindle_off(clusters{i-1}(end))) / fs;
    ITI = [ITI gap];
end

Mean_ITI = mean(ITI);

% ==============================
% OUTPUT
% ==============================
out.spindle_on = spindle_on;
out.spindle_off = spindle_off;
out.spindle_center = spindle_center;

out.frequency = sp_freq;
out.amplitude = sp_amp;

out.local_density = local_density;

out.NumCluster = NumCluster;
out.PropCluster = PropCluster;
out.AvgPerCluster = AvgPerCluster;

out.Mean_ISI = Mean_ISI;
out.Mean_ITI = Mean_ITI;

end