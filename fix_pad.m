function [new_pad] = fix_pad(pad,gap_indx)
    % considering that pad is a multiple of 4
    if mod(gap_indx-pad-1,pad) == 0
        new_pad = pad;
    else
        new_pad = pad+mod(gap_indx-pad-1,pad);
    end

    % left padding depends on first index of "gap_indx"
    if length(new_pad) >= 1
        new_pad = new_pad(1);
    end
    
end

