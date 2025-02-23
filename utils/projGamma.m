function x = projGamma(x,mask , gapped)
    %x(mask) = 0
    %mask = mask(:,1:440);
    %gapped = gapped(:,1:440);
    %gapped(mask) = 0
    %mask
    %x(mask) = gapped(mask);
    x = x.*(1-mask) + gapped.*mask;
    
end