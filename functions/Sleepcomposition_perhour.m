function [Sleep_com_perhour] = Sleepcomposition_perhour(binary_vector, start, stop, perhour)
% input: perhour -mins

    seconds_per_period = perhour * 60; % second
    periods_number = length(start:perhour:stop) - 1;


    percentage    = zeros(periods_number,1);
    bout          = zeros(periods_number,1);
    bout_duration = zeros(periods_number,1);
    [onset, offset] = binary_to_OnOff(binary_vector);
    Bout_Duration = offset - onset; 
    
    for p = 1:periods_number
        
        period_start = (p-1)*seconds_per_period + 1;
        period_end   = p*seconds_per_period;
        if period_end > length(binary_vector)
            period_end = length(binary_vector);
        end
        % percentage
        segment = binary_vector(period_start:period_end); 
        percentage(p) = sum(segment) / length(segment) * 100;

        % bout & bout duration
        idx = find( onset >= period_start & onset <= period_end );

        bout(p) = length(idx);
        bout_duration(p) = mean(Bout_Duration(idx));
    end

    Sleep_com_perhour = [percentage, bout, bout_duration];
end