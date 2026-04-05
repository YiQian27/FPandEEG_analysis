function [Positive_local_max_indices, Positive_local_max_slope] = NE_peaks_finding(data, Fs, do_plot, window_size, rebulid_level)

    % Set defaults
    if nargin < 3 || isempty(do_plot)
        do_plot = false;
    end
    if nargin < 4 || isempty(window_size)
        window_size = 4;
    end
    if nargin < 5 || isempty(rebulid_level)
        rebulid_level = 3;
    end
    
    window_samples = window_size * Fs;
    num_windows = floor((length(data) - window_samples) / Fs); 
    slope = zeros(1, num_windows);
    
    for k = 1:num_windows
        start_index = (k - 1) * Fs + 1;
        end_index = start_index + window_samples - 1;
        if end_index > length(data)
            break;
        end
        window_data = data(start_index:end_index);
        time_vector = (start_index:end_index) / Fs;
        p = polyfit(time_vector, window_data, 1); 
        slope(k) = p(1);
    end

    [c_slope, l_slope] = wavedec(slope, 30, 'db4');

    n = rebulid_level;  % Number of levels to reconstruct
    for i = 1:n
        a3(i,:) = wrcoef('a', c_slope, l_slope, 'db4', i);
    end
    slope = a3(rebulid_level,:);  % Rebuilt slope signal
    
    % Detect local maxima
    local_max_indices = islocalmax(slope, 'MinProminence', 0.05);
    local_max_slope = slope(local_max_indices);
    

    positive_mask = local_max_slope > 0;
    
    Positive_local_max_indices = find(local_max_indices);
    Positive_local_max_indices = Positive_local_max_indices(positive_mask);
    
    Positive_local_max_slope = local_max_slope(positive_mask);



    % Plot if requested
    if do_plot
        figure;
        plot(1:num_windows, slope, 'b');
        hold on
        plot(Positive_local_max_indices, Positive_local_max_slope, 'ro', 'LineWidth', 1.5);
        title('Slope Signal with Positive Local Maxima');
        xlabel('Window Index');
        ylabel('Slope');
        grid on
        hold off
    end
end