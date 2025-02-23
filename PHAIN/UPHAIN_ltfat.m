function [outsig, snr_procedure] = UPHAIN_ltfat(insig, mask, param, paramsolver, oracle)

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

rotateFlag = true;

% param.Ls = length(insig);
% param.F = frameaccel(param.F, param.Ls);

%[g,info] = gabwin({param.wtype, w},a,M,param.Ls);
%g = g/norm(g)*sqrt(a/w);
diff_win = numericalDiffWin(param.g);

%diff_win2 = numericalDiffWin(param.g);

% figure;plot(diff_win);
% figure;plot(diff_win2);
% snr(diff_win,diff_win-diff_win2)
% % DGTs (original, its adjoint, and one with the differentiated window)

% G = @(x) framecoef2tf(param.F,frana(param.F, x));
% G_adj = @(u) frsyn(param.F,frametf2coef(param.F,u));
lt = [0,1];

G = @(x) comp_sepdgtreal(x,param.g,param.a,param.M,0);
G_adj = @(u) comp_idgtreal(u,param.g,param.a,param.M,lt,0);

% DiffF = frame('dgtreal', diff_win, param.a, param.M);
% DiffF = frameaccel(DiffF, param.Ls);
% G_diff = @(x) framecoef2tf(DiffF,frana(DiffF, x));
G_diff = @(x) comp_sepdgtreal(x,diff_win,param.a,param.M,0);


omega = @(x) calcInstFreq(G(x), G_diff(x), M, w, rotateFlag);

% operator to correct phase rotation and its adjoint
R = @(z, omega) instPhaseCorrection(z, omega, a, M);
R_adj = @(z, omega) invInstPhaseCorrection(z, omega, a, M);

% time-directional difference
D = @(z) z(:,1:end-1) - z(:,2:end);
D_adj = @(z) [z(:,1), (z(:,2:end) - z(:,1:end-1)), -z(:,end)];

% iPC-DGT
hatG = @(x, omega) D(R(G(x), omega));
hatG_adj = @(u, omega) G_adj(R_adj(D_adj(u), omega));


%%

soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);
param.proj = @(x) projGamma(x, mask, insig.*mask);
param.G = @(x) G(x);

hold on;
if strcmp(param.type,'U')
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    param.prox = @(z) soft(z, lambda/sigma);

    omega_y = omega(insig);
    param.L = @(x) hatG(x, omega_y);
    param.L_adj = @(u) hatG_adj(u, omega_y);

    paramsolver.x0 = insig;
    paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
    x_old = insig;

    snr_procedure = NaN(paramsolver.I, paramsolver.J);

    for j = 1:paramsolver.J

        [x_hat, snr_procedure(:, j)] = CP(param, paramsolver, oracle, mask);
        
        if norm(x_old - x_hat) < paramsolver.epsilon
            break
        end

        omega_x_hat = omega(x_hat);
        param.L = @(x) hatG(x, omega_x_hat);
        param.L_adj = @(u) hatG_adj(u, omega_x_hat);

        x_old = x_hat;

    end
    outsig = x_hat;
end
