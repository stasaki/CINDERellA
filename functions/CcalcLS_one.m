

function [localscore] = CcalcLS_one(expdata,gene_indx,pasets_cand)
nSamples = size(expdata,2);
y = expdata(gene_indx,:)';

% calc local score
if isempty(pasets_cand)
    
    X1=ones(nSamples,1);
    beta = X1 \ y;
    
    yhat = X1 *beta;
    residuals = y - yhat;
    params(2)=sqrt(residuals'*residuals/nSamples);
    
    loglik = -normlike(params,residuals);
    localscore = loglik - 0.5*length(beta)*log(nSamples);
else
    X = [ones(1,nSamples);expdata(pasets_cand,:)]';
    
    beta = X \ y;
    
    yhat = X*beta;
    residuals = y - yhat;
    loglik = -my_normlike(residuals,0,sqrt(residuals'*residuals/nSamples));
    localscore = loglik - 0.5*length(beta)*log(nSamples);% should be replaced by my_normlike
end

