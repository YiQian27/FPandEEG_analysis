function [infra_fre, spectrum_epoch_weighted] = analyze_ISO_spectrum(signal_trace, fs_ref, NREMinclMA_periods, min_period_dur, wavelet_number, do_plot)
% ANALYZE_ISO_SPECTRUM Analyze infraslow oscillations in sigma power using Morlet wavelet transform
%
% Inputs:
%   signal_trace - vector of the signal (e.g., sigma power time series)
%   fs_ref - sampling frequency
%   NREMinclMA_periods - Nx2 matrix of start/end times (in seconds)
%   min_period_dur - (optional) minimum duration of each valid period (s) [default: 240]
%   wavelet_number - (optional) wavelet width parameter [default: 5]
%   do_plot - (optional) whether to plot the result [default: true]
%
% Outputs:
%   infra_fre - frequency vector (Hz)
%   spectrum_epoch_weighted - weighted average power spectral density across valid periods

    % Set defaults
    if nargin < 4 || isempty(min_period_dur)
        min_period_dur = 120;
    end
    if nargin < 5 || isempty(wavelet_number)
        wavelet_number = 5;
    end
    if nargin < 6
        do_plot = true;
    end

    t1 = NREMinclMA_periods(:,1) + 1;
    t2 = NREMinclMA_periods(:,2);

    tsamp1 = max(floor(t1 * fs_ref) + 1, 1);
    tsamp2 = min(floor(t2 * fs_ref), length(signal_trace));

    infra_fre = 0.004:0.002:0.12;
    spectrum_epoch = [];
    period_duration = [];

    for i = 1:numel(tsamp1)
        period_length_i = tsamp2(i) - tsamp1(i);
        if period_length_i < min_period_dur * fs_ref
            continue
        end

        tsamp2(i) = min(tsamp2(i), length(signal_trace));
        period_duration = [period_duration, period_length_i / fs_ref];

        NREM_data = signal_trace(tsamp1(i):tsamp2(i));
        [p, s, mu] = polyfit((1:numel(NREM_data))', NREM_data, 5);
        f_y = polyval(p, (1:numel(NREM_data))', [], mu);
        detrend_data = NREM_data - f_y';

        conv_result_fft = ISr_morlet_transform( ...
            detrend_data, infra_fre, wavelet_number, -100:1/fs_ref:100, 0, fs_ref);

        spectrum_epoch = [spectrum_epoch, mean(abs(conv_result_fft).^2, 2)];
    end

    spectrum_epoch_weighted = sum(period_duration .* spectrum_epoch, 2) / sum(period_duration);

    if do_plot
        figure('Position', [100 100 800 400])
        plot(infra_fre, spectrum_epoch_weighted, 'LineWidth', 1.5)
        xlabel('Frequency (Hz)', 'FontSize', 12)
        ylabel('Power Spectral Density', 'FontSize', 12)
        title('Infraslow Oscillations in Sigma Power', 'FontSize', 14)
        grid on
    end

end