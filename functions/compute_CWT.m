function [power_tf, f, t_eeg_ds] = compute_CWT(EEG, f_target, EEG_fs, binary_vector)
% do-sample to 100 Hz
fs_ds = 100;  
t_original = (0:length(EEG)-1)/EEG_fs;
t_ds = 0:1/fs_ds:t_original(end);
EEG_ds = interp1(t_original, EEG, t_ds, 'linear');
t_eeg_ds = (0:length(EEG_ds)-1)/fs_ds;
% NREM mask, only in NREM
NREM_mask_ds = repelem(binary_vector, fs_ds);
NREM_mask_ds = NREM_mask_ds(1:length(EEG_ds));
NREM_mask_ds = logical(NREM_mask_ds);

[cfs, f] = cwt(EEG_ds, 'amor', fs_ds, ...
               'FrequencyLimits', [min(f_target) max(f_target)], ...
               'VoicesPerOctave', 10);

f = f(:); 
idx_f = interp1(f, 1:length(f), f_target, 'nearest', 'extrap');
cfs = cfs(idx_f, :);
f = f_target;

power_tf = abs(cfs).^2;  % [freq x time]

power_tf(:,~NREM_mask_ds) = NaN;

end