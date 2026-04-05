function [t_slope, slope] = derivative (signal_trace, Fs, window_size, step_size, plot_command)

% window size: s
% step size: s

Fs = Fs; 
data = signal_trace; 
T = length(data)./Fs ; 

window_samples = window_size * Fs;
step_samples = step_size * Fs;

num_windows = floor((length(data) - window_samples) / step_samples); 

slope = zeros(1, num_windows);

for k = 1:num_windows
    start_index = round((k - 1) * step_samples) + 1;
    end_index = round(start_index + window_samples) - 1;
    
    if end_index > length(data)
        break;
    end
    
    window_data = data(start_index:end_index);
    time_vector = (start_index:end_index) / Fs;
    p = polyfit(time_vector, window_data, 1); 
    slope(k) = p(1);
end

% 
% [c_slope,l_slope]=wavedec(slope,30,'db4');
% 
% % figure;
% n=4;
% clear a3;
% clear b3;
% for i=1:n
%     a3(i,:)=wrcoef('a',c_slope,l_slope,'db4',i);
%     b3(i,:)=wrcoef('d',c_slope,l_slope,'db4',i);
%     % subplot(n,2,2*i-1);plot(t,a3(i,:));
%     % subplot(n,2,2*i);plot(t,b3(i,:));
% end
% 
% y1=a3(1,:);
% 
if plot_command == true
    figure;
    plot(slope);
    title('slope');

    % subplot(2,1,2);
    % plot(y1);
    % title('Denoised slope');
end 
% 
% 
% slope = y1;

t_slope = linspace (1,T, length(slope)); 


