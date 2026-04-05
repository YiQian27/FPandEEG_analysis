
function index = get_cor_indices(period1, period2)
% 
%
% input:
%   offset - the one period  
%   onset - the other period
%
% output: index of the corresponding

    index = []; 
    for i = 1:size(period1,1)
        start = period1(i,2);
        idx = find(period2(:,1) == start);
        if ~isempty(idx)
            index(end + 1,:) = [i,idx];
        end
    end
end