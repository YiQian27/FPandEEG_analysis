function [new_time,data_out] = move_prctile(t,y,buffer_size,overlap,p)
% NB : buffer size and overlap are integer numbers (samples)
% data (in , out) are 1D arrays (vectors)
shift = buffer_size-overlap;    % nb of samples between 2 contiguous buffers  
samples = numel(y);
nb_of_loops = fix((samples-buffer_size)/shift +1);
    for k=1:nb_of_loops
        start_index = 1+(k-1)*shift;
        stop_index = min(start_index+ buffer_size-1,samples);
        x_index(k) = round((start_index+stop_index)/2);  
        data_out(k) = prctile(y(start_index:stop_index),p);  % 
    end
new_time = t(x_index); % time values are computed at the center of the buffer
end

