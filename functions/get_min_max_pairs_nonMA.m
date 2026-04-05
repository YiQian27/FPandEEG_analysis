function Range_pair_sec_nonMA = get_min_max_pairs_nonMA( ...
    min_locs, min_values, ...
    peak_locs, peaks, ...
    NREMinclMA_periods, MA_periods, ...
    ds_signal_fs, ds_delta465_filt_1)

% -----------------------------
% Step 1: Get min/max in NREM (incl MA)
% -----------------------------
locs_min_NREMinclMA = cell2mat(arrayfun(@(i) ...
    min_locs(min_locs >= (NREMinclMA_periods(i,1)+1)*ds_signal_fs & ...
             min_locs <  (NREMinclMA_periods(i,2)  )*ds_signal_fs), ...
    1:size(NREMinclMA_periods,1), 'UniformOutput', false));

locs_max_NREMinclMA = cell2mat(arrayfun(@(i) ...
    peak_locs(peak_locs >= (NREMinclMA_periods(i,1)+1)*ds_signal_fs & ...
              peak_locs <  (NREMinclMA_periods(i,2)  )*ds_signal_fs), ...
    1:size(NREMinclMA_periods,1), 'UniformOutput', false));

locs_min = sort(locs_min_NREMinclMA);
locs_max = sort(locs_max_NREMinclMA);

vals_min = min_values(ismember(min_locs, locs_min));
vals_max = peaks(ismember(peak_locs, locs_max));

% -----------------------------
% Step 2: Pair min → max
% -----------------------------
window_sec = 20;
window_samples = window_sec * ds_signal_fs;

tmp_pairs = [];

for i = 1:length(locs_max)

    max_loc = locs_max(i);
    max_val = vals_max(i);

    if i == 1
        lower_bound = max_loc - window_samples;
    else
        lower_bound = max(locs_max(i-1), max_loc - window_samples);
    end

    idx_min_valid = find( ...
        locs_min < max_loc & ...
        locs_min >= lower_bound);

    if isempty(idx_min_valid)
        continue
    end

    [min_val, idx_low] = min(vals_min(idx_min_valid));
    min_loc = locs_min(idx_min_valid(idx_low));

    tmp_pairs(end+1,:) = [min_loc, max_loc, min_val, max_val];
end

% -----------------------------
% Step 3: Unique min handling
% -----------------------------
final_pairs = [];
unique_mins = unique(tmp_pairs(:,1));

for i = 1:length(unique_mins)

    min_loc = unique_mins(i);
    idx = tmp_pairs(:,1) == min_loc;

    if sum(idx) == 1
        final_pairs(end+1,:) = tmp_pairs(idx,:);
    else
        [~, idx_maxbest] = max(tmp_pairs(idx,4));
        rows = find(idx);
        final_pairs(end+1,:) = tmp_pairs(rows(idx_maxbest),:);
    end
end

pair_min_locs = final_pairs(:,1);
pair_max_locs = final_pairs(:,2);

% -----------------------------
% Step 4: Refine on raw signal
% -----------------------------
fs = ds_signal_fs;
search_win = round(1 * fs);

raw_signal = ds_delta465_filt_1;
N = length(raw_signal);

raw_min_locs = nan(size(pair_min_locs));
raw_max_locs = nan(size(pair_max_locs));

for i = 1:length(pair_min_locs)

    % MIN
    idx = max(1, pair_min_locs(i)-search_win) : min(N, pair_min_locs(i)+search_win);
    [~, idx_rel] = min(raw_signal(idx));
    raw_min_locs(i) = idx(idx_rel);

    % MAX
    idx = max(1, pair_max_locs(i)-search_win) : min(N, pair_max_locs(i)+search_win);
    [~, idx_rel] = max(raw_signal(idx));
    raw_max_locs(i) = idx(idx_rel);
end

assert(all(raw_max_locs > raw_min_locs), 'Max must be after Min');

% -----------------------------
% Step 5: Convert to seconds
% -----------------------------
Range_pair_loc = [raw_min_locs, raw_max_locs];
Range_pair_sec = Range_pair_loc ./ ds_signal_fs;

% -----------------------------
% Step 6: Remove overlap with MA
% -----------------------------
keep_idx = true(size(Range_pair_sec,1),1);

for i = 1:size(Range_pair_sec,1)

    pair_start = Range_pair_sec(i,1);
    pair_end   = Range_pair_sec(i,2);

    overlap = any( ...
        (pair_start <= MA_periods(:,2)) & ...
        (pair_end   >= MA_periods(:,1)) );

    if overlap
        keep_idx(i) = false;
    end
end

Range_pair_sec_nonMA = Range_pair_sec(keep_idx,:);

end