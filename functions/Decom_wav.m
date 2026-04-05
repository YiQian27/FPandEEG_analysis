function [signal_clean] = Decom_wav( signal, use_lvl, do_plot)

% INPUT:
%   use_lvl            - level
%   do_plot            - plot or not
%
    ds_sec_signal = 1:1:length(signal);
    start_idx = 11;
    wavelet  = 'db4';
    nlevel   = 30;
    s = signal(start_idx:end);
    t_out = ds_sec_signal(start_idx:end);
    y = wden(s, 'heursure', 's', 'sln', 5, wavelet);
    [c, l] = wavedec(y, nlevel, wavelet);
    y1 = wrcoef('a', c, l, wavelet, use_lvl);

    signal_clean = [signal(1:10), y1];

    %% Plot
    if do_plot
        figure;
        plot(t_out, y1, 'LineWidth', 1.2);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title(['Wavelet approximation level a', num2str(use_lvl)]);
        box off;
    end

end