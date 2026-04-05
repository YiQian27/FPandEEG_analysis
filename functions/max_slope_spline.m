function [x_at_max_slope, y_at_max_slope, max_slope, f] = max_slope_spline(t, delta_trace, t_start_idx, t_end_idx, smoothing_param)
%MAX_SLOPE_SPLINE Fit a smoothing spline and find the max slope in a given interval
%
% Inputs:
%   t             - time vector
%   delta_trace   - data vector corresponding to t
%   t_start_idx   - starting index in t for region of interest
%   t_end_idx     - ending index in t for region of interest
%   smoothing_param - smoothing parameter for the spline (0~1)
%
% Outputs:
%   x_at_max_slope - time at max slope within ROI
%   y_at_max_slope - spline value at max slope
%   max_slope      - maximum slope value
%   f              - fitted spline object

    % Fit smoothing spline
    f = fit(t(:), delta_trace(:), 'smoothingspline', 'SmoothingParam', smoothing_param);

    % Evaluate spline and derivative
    xx = linspace(min(t), max(t), 2000)';
    yy = feval(f, xx);
    dy = differentiate(f, xx);

    % Define region of interest
    t_start = t(t_start_idx);
    t_end   = t(t_end_idx);
    roi_idx = find(xx >= t_start & xx <= t_end);

    % Find max slope in ROI
    dy_roi = dy(roi_idx);
    [max_slope, local_max_idx] = max(dy_roi);
    x_at_max_slope = xx(roi_idx(local_max_idx));
    y_at_max_slope = yy(roi_idx(local_max_idx));

    % Plot
    figure; hold on;
    plot(t, delta_trace, 'k.', 'DisplayName', 'Raw data');
    plot(xx, yy, 'b-', 'LineWidth', 2, 'DisplayName', 'Spline fit');
    plot(x_at_max_slope, y_at_max_slope, 'ro', 'MarkerSize', 8, ...
        'LineWidth', 2, 'DisplayName', 'Max slope point in ROI');
    xlabel('Time');
    ylabel('Delta Trace');
    title('Spline Fit with Maximum Slope in Specified Interval');
    legend('Location', 'best');
    grid on;
end