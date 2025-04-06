function [outsig] = PHAINmain_TF(insig, mask, param, paramsolver,oracle)

% param
%   .a
%   .M
%   .w
%   .type 


% paramsolver
%   .sigma .... step size
%   .tau ...... step size
%   .alpha .... relaxation paramter
%   .lambda ... threshold
%   .epsilon .. stop criterion
%   .x0 ....... initial value of primal variable
%   .u0 ....... initial value of dual variable
%   .I ........ number of inner iterations
%   .J ........ number of outer iterations


%% iPC DGT

a = param.a;
M = param.M;
w = param.w;

[win, ~] = generalizedCosWin(w, 'hanning');
tight_win = calcCanonicalTightWindow(win, a);
tight_win = tight_win/norm(tight_win)*sqrt(a/w);
diff_win = numericalDiffWin(tight_win);
    
zeroPhaseFlag = true;
rotateFlag = true;

signal_l = size(insig,2)* a;
[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(signal_l, w, a, M);

% DGTs (original, its adjoint, and one with the differentiated window)
G = @(x) FDGT(x, tight_win, sigIdx, M, rotIdx, zeroPhaseFlag);
G_adj = @(u) invFDGT(u, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*w;
G_diff = @(x) FDGT(x, diff_win, sigIdx, M, rotIdx, zeroPhaseFlag);

% function to calculate the instantaneous frequency of the input signal
omega = @(x) calcInstFreq(x, G_diff(G_adj(x)), M, w, rotateFlag);

% operator to correct phase rotation and its adjoint
R = @(z, omega) instPhaseCorrection(z, omega, a, M);
R_adj = @(z, omega) invInstPhaseCorrection(z, omega, a, M);

% time-directional difference (time variation)
D = @(z) z(:,1:end-1) - z(:,2:end);
D_adj = @(z) [z(:,1), (z(:,2:end) - z(:,1:end-1)), -z(:,end)];

param.G = G;
param.G_adj = G_adj;
%%

soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);
param.proj = @(x) projGamma(x, mask, insig.*mask);

if strcmp(param.type,'U')

    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    param.prox = @(z) soft(z, lambda/sigma);

    % first calculate omega from the input signal
    omega_y = omega(insig);
    param.L = @(x) D(R(x, omega_y));
    param.L_adj =@(u) R_adj(D_adj(u), omega_y);
    x_old = insig;
    
    % set starting x and u for CP alg.
    paramsolver.x0 = insig;
    paramsolver.u0 = zeros(size(param.L(insig)));

    for j = 1:paramsolver.J
        % calculate CP
        x_hat = CP_TF(param, paramsolver,oracle);

        % stopping criterion
        if norm(x_old - x_hat) < paramsolver.epsilon 
            break 
        end

        % calculate new IF and redefine L and L_adj
        omega_x_hat = omega(x_hat);
        param.L = @(x) D(R(x, omega_x_hat));
        param.L_adj =@(u) R_adj(D_adj(u), omega_x_hat);
        x_old = x_hat;
    end
    outsig = x_hat;
end
