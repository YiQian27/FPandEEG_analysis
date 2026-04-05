function out = iff(cond, val_true, val_false)
%IFF  Inline if-else function
%   out = iff(cond, a, b)
%   if cond is true, out = a; otherwise out = b

    if cond
        out = val_true;
    else
        out = val_false;
    end
end