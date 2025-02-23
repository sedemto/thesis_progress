%function [x] = CP_TF(param, paramsolver, oracle, mask)
function [x] = CP_TF(param, paramsolver)

%% initialization

x = paramsolver.x0;
u = paramsolver.u0;

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;


%% iteration

for i = 1:paramsolver.I

    p = param.proj(x - tau*sigma*param.L2_adj(u));
    v = u+param.L2(2*p - x);
    q = v - param.prox(v);
    x = x + alpha*(p - x);
    u = u + alpha*(q - u);
    
end
