function x = projGamma(x,mask , gapped)
    x = x.*(1-mask) + gapped.*mask;
    
end