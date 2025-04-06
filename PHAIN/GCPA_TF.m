%function [x] = CP_TF(param, paramsolver, oracle, mask)
function [x] = GCPA_TF(param, paramsolver, oracle)

%% initialization

x = paramsolver.x0;
u = paramsolver.u0;
v = paramsolver.v0;

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;
eta = paramsolver.eta;


SNRS = zeros(4,100);

%% iteration

for i = 1:paramsolver.I

    tmp_r = v + eta*param.G(x-tau*(param.L_adj(u)+param.G_adj(v)));
    r = tmp_r - eta*param.proj(tmp_r/eta);
    
    p = x - tau*(param.L_adj(u)+param.G_adj(r));

    tmp_q = u + sigma*param.L(2*p-x);
    q = tmp_q - sigma*param.prox(tmp_q/sigma);

%     tmp_r = v + eta*param.L(x-tau*(param.G_adj(u)+param.L_adj(v)));
%     r = tmp_r - eta*param.prox(tmp_r/eta);
%     
%     p = x - tau*(param.G_adj(u)+param.L_adj(r));
% 
%     tmp_q = u + sigma*param.G(2*p-x);
%     q = tmp_q - sigma*param.proj(tmp_q/sigma);


%     v = u+sigma*param.L(x);
%     q = v - sigma*param.proj(v/sigma);
% 
%     p = param.prox(x-tau*param.L_adj(2*q-u));
    
    x = x + alpha*(p - x);
    u = u + alpha*(q - u);
    v = v + alpha*(r - v);

    
%     x = param.L_adj(param.L(x));
    x_spec = param.proj(param.G(x));
    x = param.G_adj(x_spec);

%     objective = param.f(param.G(x))+1000000*param.g(x);%+norm(param.proj(param.L(x)));
%     snr_iter = snr(param.G_adj(oracle), param.G_adj(oracle)-param.G_adj(x_spec));
%     snr_gap  = snr(oracle(:,param.pad+1),oracle(:,param.pad+1)-x_spec(:,param.pad+1));
%     snr_iter2 = snr(oracle, oracle-x_spec);
%     SNRS(1,i) = snr_iter;
%     SNRS(2,i) = snr_iter2;
%     SNRS(3,i) = snr_gap;
%     SNRS(4,i) = objective;
    
    
end
% figure, imagesc(20*log10(abs(oracle))), colorbar, axis xy
% figure, imagesc(angle(oracle)), colorbar, axis xy
% 
% figure, imagesc(20*log10(abs(x_spec))), colorbar, axis xy
% figure, imagesc(angle(x_spec)), colorbar, axis xy
% figure; plot(param.G_adj(oracle),'Color',[0 0 1 1]), hold on, plot(param.G_adj(x_spec),'Color',[1 0 0 0.5])
% disp("here1")