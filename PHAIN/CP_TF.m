%function [x] = CP_TF(param, paramsolver, oracle, mask)
function [x] = CP_TF(param, paramsolver, oracle)

%% initialization

x = paramsolver.x0;
u = paramsolver.u0;

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;

SNRS = zeros(2,100);

%% iteration

for i = 1:paramsolver.I
%     x = param.G(param.G_adj(x));
    p = param.proj(x - tau*sigma*param.L_adj(u));
    v = u+param.L(2*p - x);
    q = v - param.prox(v);
    x = x + alpha*(p - x);
    u = u + alpha*(q - u);
    snr_iter = snr(oracle, oracle-x);
    snr_gap  = snr(oracle(:,8+1),oracle(:,8+1)-x(:,8+1));
    SNRS(1,i) = snr_iter;
    SNRS(2,i) = snr_gap;
end
disp("here1")