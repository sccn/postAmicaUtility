function LLt = getmodLLt(x,modl,N)
% function LLt = getmodLLt(x,mod)
% Get the (multiple model) likelihood of data x under the model(s) in the
% AMICA output mod
c = modl.c;
W = modl.W;
S = modl.S(1:modl.num_pcs,:);
gm = modl.mod_prob;
alpha = modl.alpha;
mu = modl.mu;
sbeta = modl.sbeta;
rho = modl.rho;
mn = modl.data_mean;

[nx,T] = size(x);
n = size(c,1);
M = length(gm);
m = size(alpha,1);
if nargin < 3
    N = 10000;
end

minlog = -1500;

num_blocks = floor(T/N);
lastblocksize = N + mod(T,N);

N2 = lastblocksize;

LLt = zeros(M,T);

btmp = zeros(n,N2);
ytmp = zeros(1,N2);

strt = 1;

% sphere the data
for i = 1:nx
    x(i,:) = x(i,:) - mn(i);
end

for k = 1:num_blocks
    xstrt = (k-1)*N + 1;
    if k < num_blocks
        xstp = k*N;
    else
        xstp = T;
    end
    x(1:n,xstrt:xstp) = S * x(:,xstrt:xstp);
end

for h = 1:M
    LLt(h,:) = 0.5*log(abs(det(S*S'))) + log(abs(det(W(:,:,h)))) + log(gm(h));
end


for k = 1:num_blocks
    if mod(k,1)==0
        disp(['doing block ' int2str(k) ' ...']); pause(0.2);
    end
    xstrt = (k-1)*N + 1;
    if k < num_blocks
        xstp = k*N;
        stp = N;
    else
        xstp = T;
        stp = lastblocksize;
    end

    for h = 1:M
        btmp(:,strt:stp) = W(:,:,h) * x(1:n,xstrt:xstp);
        for i = 1:n
            for j = 1:m
                ytmp(strt:stp) = sbeta(j,i,h) * ( btmp(i,strt:stp) - mu(j,i,h) - c(i,h));
                %Q(j,strt:stp) = log(alpha(j,i,h)) + log(sbeta(j,i,h)) + feval('logpfun',ytmp(strt:stp),rho(j,i,h));
                if alpha(j,i,h) > 0
                    Q(j,strt:stp) = log(alpha(j,i,h)) + log(sbeta(j,i,h)) - abs(ytmp(strt:stp)).^rho(j,i,h) - log(2) - gammaln(1+1/rho(j,i,h));               
                else
                    Q(j,strt:stp) = minlog;
                end
            end
            Qmax(:,strt:stp) = ones(m,1) * max(Q(:,strt:stp),[],1);
            LLt(h,xstrt:xstp) = LLt(h,xstrt:xstp) + Qmax(1,strt:stp) + log(sum(exp(Q(:,strt:stp) - Qmax(:,strt:stp)),1));
        end    
    end
end



function f = pfun(x,rho)
f = (1 / (2*gamma(1+1/rho))) * exp( - abs(x).^rho );

function f = logpfun(x,rho)
f = -abs(x).^rho - log(2) - gammaln(1+1/rho);
