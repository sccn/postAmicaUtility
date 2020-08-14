function [MI,T,Hu,hu,Hx,hx] = get_mi(u,x,nlags,p,nbins,verb)
% function [MI,T,Hu,hu,Hx,hx] = get_mi(u,x,nlags,p,nbins,verb)
%
% Calculate the pairwise mutual information between rows of u, or between
% rows of u and rows of x if second argument is a matrix. If nlags > 0,
% calculate for lags -nlags:nlags. The output is a matrix of size
% [nu,nx,2*nlags+1] of pairwise mutual information values, with zero lag
% at MI(:,:,nlags+1). Uses nbins bins for pdf of u(i,:), equally spaced
% between mean +/- 5 std dev. (or data min/max), similarly for pdfs of x rows.
%
% Inputs:
%           u       Matrix (nu by Nu) of nu time series.
%
%           x       Optional second input matrix (nx by Nx) of nx time
%                   series. Nx must equal Nu, i.e. all time series must
%                   have the same length. Enter 0 to bypass in order to
%                   set nlags or nbins.
%
%           nlags   Number of lags (in addition to zero lag) to compute
%                   mutual information (before and after zero) lag.
%                   Default is 0.
%
%           p       p value (right tail probability) for determining the
%                   threshold T that is returned for rejecting hypothesis
%                   of independence (MI==0).
%
%           nbins   Number of bins to use in computing pdfs. Default is
%                   round(3*log2(1+Nu/10)). Enter 0 to use default.
%
%           verb    Verbose flag. 1 = verbose (default), 0 = no text output.
%
% Outputs:
%           MI      If no second input x, or x==0, then MI is an
%                   (nu by nu by 2*nlags+1) matrix of pairwise mutual
%                   information between rows of u at lags -nlags to nlags.
%                   If both inputs u and x are present, then MI is an
%                   (nu by nx by 2*nlags+1) matrix of pairwise mutual
%                   information between rows of u and rows of x. For
%                   example, MI(i,j,nlags+1+g) is the mutual information
%                   between u(i,1+g:N+g) and u(j,1:N), or between
%                   u(i,1+g:N+g) and x(j,1:N). So a peak at g+nlags+1 in
%                   MI(i,j,:) means i lags j by g time points (j leads i by
%                   g time points.) Diagonal entries are zero.
%
%           T       If significance input p is given, then T is threshold at
%                   which hypothesis of independence is rejected at level p
%                   signifcance. If no in is input, or p=0, then T is array
%                   of significance thresholds for p = [0.001 0.01 0.05 0.1 0.5].
%                   Significance is based treating MI as a chi squared random
%                   variable. A threshold can be determined at the "p"
%                   significance level using:
%
%                       >> T = chi2inv(1-p,dof) / (2*Nu);
%
%                   By the Likelihood Ratio Test theorem, 2*Nu*MI is chi
%                   squared with nbins*(nbins-2) d.o.f. So MI has a mean of
%                   nbins*(nbins-2)/(2*Nu) for independent r.v.'s,
%                   If you don't have access to chi2inv(), a simple threshold
%                   of dof/Nu can be used.
%
%           Hu      Vector nu by 1 of discrete histogram entropies of rows of u.
%           hu      Vector nu by 1 of differential entropies of rows of u.
%
%           Hx      Vector nu by 1 of discrete histogram entropies of rows of x.    
%           hx      Vector nu by 1 of differential entropies of rows of x.           
%
varrng = 1;
nstddev = 5;
if nargin < 6
    verb = 1;
end
[nu,Nu] = size(u);
if (nargin < 5) || (isempty(nbins)) || (nbins == 0)
    nbins = round(3*log2(1+Nu/10));
end
if (nargin < 3) || (nlags == 0) || (isempty(nlags))
    nlags = 0;
end
if nargin >= 2
    [nx,Nx] = size(x);
    if Nx == 1
        MI = zeros(nu,nu,2*nlags+1); 
        dox = 0;
    else        
        if Nx ~= Nu
            error('input time series must be same length');
        end
        MI = zeros(nu,nx,2*nlags+1);
        dox = 1;
    end
else
    MI = zeros(nu,nu,2*nlags+1);
    dox = 0;
end

% bin the time series
if verb
    disp('binning the time series ...'); pause(0.1);
end
Hu = zeros(nu,2*nlags+1);
deltau = zeros(nu,1);
for i = 1:nu
    if varrng
        um = mean(u(i,:));
        us = std(u(i,:));
        umax = min(max(u(i,:)),um + nstddev*us);
        umin = max(min(u(i,:)),um - nstddev*us);
    else
        umax = max(u(i,:));
        umin = min(u(i,:));
    end
    deltau(i) = (umax-umin)/nbins;
    u(i,:) = 1 + round((nbins - 1) * (u(i,:) - umin) / (umax - umin));
    u(i,:) = min(nbins,u(i,:));
    u(i,:) = max(1,u(i,:));
    
    if nlags == 0
        pmfr = diff([0 find(diff(sort(u(i,:)))) Nu])/Nu;
        Hu(i) = -sum(pmfr.*log(pmfr));
    else
        for g = 0:nlags
            pmfr = diff([0 find(diff(sort(u(i,1+g:Nu)))) (Nu-g)])/(Nu-g);    
            Hu(i,nlags+1+g) = -sum(pmfr.*log(pmfr));
            if g > 0
                pmfr = diff([0 find(diff(sort(u(i,1:Nu-g)))) (Nu-g)])/(Nu-g);    
                Hu(i,nlags+1-g) = -sum(pmfr.*log(pmfr));
            end
        end
    end
end
if dox == 1
    Hx = zeros(nx,2*nlags+1);
    deltax = zeros(nx,1);
    for i = 1:nx
        if varrng
            xm = mean(x(i,:));
            xs = std(x(i,:));
            xmax = min(max(x(i,:)),xm + nstddev*xs);
            xmin = max(min(x(i,:)),xm - nstddev*xs);
        else
            xmax = max(x(i,:));
            xmin = min(x(i,:));
        end
        deltax(i) = (xmax-xmin)/nbins;
        x(i,:) = 1 + round((nbins - 1) * (x(i,:) - xmin) / (xmax - xmin));

        if nlags == 0
            pmfr = diff([0 find(diff(sort(x(i,:)))) Nx])/Nx;                        
            Hx(i) = -sum(pmfr.*log(pmfr));
        else
            for g = 0:nlags
                pmfr = diff([0 find(diff(sort(x(i,1+g:Nx)))) (Nx-g)])/(Nx-g);
                Hx(i,nlags+1+g) = -sum(pmfr.*log(pmfr));
                if g > 0
                    pmfr = diff([0 find(diff(sort(x(i,1:Nx-g)))) (Nx-g)])/(Nx-g);    
                    Hx(i,nlags+1-g) = -sum(pmfr.*log(pmfr));
                end
            end
        end
    end
else
    Hx = 0;
end

% get pairwise histograms at each lag
if verb
    disp('getting pairwise mutual information ...'); pause(0.1);
end
if dox == 0
    for i = 1:nu
        if verb
            if i == 1
                disp(['doing ' int2str(i) ' ...']); pause(0.01);
                tic;
            else
                t1 = toc;
                if nlags == 0
                    tpe = t1 / (nu-i+1);
                    trem = tpe * (nu-i)*(nu-i-1)/2;
                else
                    tpe = t1 / (nu-i+2);
                    trem = tpe * (nu-i+1)*(nu-i)/2;
                end
                disp(['doing ' int2str(i) ' ... time rem: ' num2str(trem/60) ' m']); pause(0.01);
                tic;
            end
        end
        for j = (i+1):nu
            if nlags == 0
                pmf2r = diff([0 find(diff(sort(u(i,:)+nbins*(u(j,:)-1)))) Nu])/Nu;
                H2 = -sum(pmf2r.*log(pmf2r));
                MI(i,j) = Hu(i) + Hu(j) - H2;
                MI(j,i) = MI(i,j);
            else
                % do zero lag
                pmf2r = diff([0 find(diff(sort(u(i,1:Nu)+nbins*(u(j,1:Nu)-1)))) Nu])/Nu;
                H2 = -sum(pmf2r.*log(pmf2r));
                MI(i,j,nlags+1) = Hu(i,nlags+1) + Hu(j,nlags+1) - H2;
                MI(j,i,nlags+1) = MI(i,j,nlags+1);

                % do +g and -g lags
                for g = 1:nlags
                    pmf2r = diff([0 find(diff(sort(u(i,1+g:Nu)+nbins*(u(j,1:(Nu-g))-1)))) Nu-g])/(Nu-g);
                    H2 = -sum(pmf2r.*log(pmf2r));
                    MI(i,j,nlags+1+g) = Hu(i,nlags+1+g) + Hu(j,nlags+1-g) - H2;                    
                    MI(j,i,nlags+1-g) = MI(i,j,nlags+1+g);

                    pmf2r = diff([0 find(diff(sort(u(j,1+g:Nu)+nbins*(u(i,1:(Nu-g))-1)))) Nu-g])/(Nu-g);
                    H2 = -sum(pmf2r.*log(pmf2r));
                    MI(i,j,nlags+1-g) = Hu(i,nlags+1+g) + Hu(j,nlags+1-g) - H2;                    
                    MI(j,i,nlags+1+g) = MI(i,j,nlags+1-g);
                end
            end        
        end
    end    
else
    for i = 1:nu
        if verb
            if i == 1                
                disp(['doing ' int2str(i) ' ...']); pause(0.01);
                tic;
            else
                t1 = toc; t0 = t1/nx; trem = t0 * nx * (nu-i+1);
                disp(['doing ' int2str(i) ' ... time rem: ' num2str(trem/60) ' m']); pause(0.01);
                tic;
            end
        end
        for j = 1:nx
            if nlags == 0
                pmf2r = diff([0 find(diff(sort(u(i,:)+nbins*(x(j,:)-1)))) Nu])/Nu;            
                H2 = -sum(pmf2r.*log(pmf2r));
                MI(i,j) = Hu(i) + Hx(j) - H2;
            else
                for g = 0:nlags
                    pmf2r = diff([0 find(diff(sort(u(i,1+g:Nu)+nbins*(x(j,1:(Nu-g))-1)))) Nu-g])/(Nu-g);     
                    H2 = -sum(pmf2r.*log(pmf2r));
                    MI(i,j,nlags+1+g) = Hu(i,nlags+1+g) + Hx(j,nlags+1-g) - H2;
                
                    if g > 0                
                        pmf2r = diff([0 find(diff(sort(u(i,1:(Nu-g))+nbins*(x(j,1+g:Nu)-1)))) Nu-g])/(Nu-g);        
                        H2 = -sum(pmf2r.*log(pmf2r));
                        MI(i,j,nlags+1-g) = Hu(i,nlags+1-g) + Hx(j,nlags+1+g) - H2;
                    end
                end
            end        
        end        
    end
end
Hu = Hu(:,nlags+1);
hu = Hu;
for i = 1:nu
    hu(i) = hu(i) + log(deltau(i)); 
end
if dox == 1
    Hx = Hx(:,nlags+1);
    hx = Hx;
    for i = 1:nx
        hx(i) = hx(i) + log(deltax(i));
    end
end
dof = nbins*(nbins-2);
if ~isempty(ver('stats'))
    if (nargin >= 4) && (p > 0)
        T = chi2inv(1-p,dof)/(2*Nu);
    else
        T = chi2inv([0.999 0.99 0.95 0.9 0.5],dof)/(2*Nu);
    end
else
    T = dof/Nu;
    if verb
        disp('no stats toolbox detected');
    end
end
if verb
    disp(['Degrees of freedom = ' int2str(dof)]);
    disp(['Thresholds = ' num2str(T) ' (p=0.001, 0.01, 0.05, 0.1, 0.5)']);
end