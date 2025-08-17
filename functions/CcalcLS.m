function [LS] = CcalcLS(expdata, varargin)
% CcalcLS pre-computing local score for network fragments
% INPUTS:
%	expdata:	gene expression data
%   pa_limit: the number of parents allowed
%	hprior:     matrix of acceptable edges: 1-allowed, 0-banned
%   sprior:   log of prior for the egdge
% OUTPUTS:
%	LS:		local score object
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License



% rng('shuffle')
[nGenes, nSamples] = size(expdata);

% set default paramters
pa_limit=3;
Sprior = zeros(nGenes);
Hprior=[];

args = varargin;
nargs = length(args);

% set paramters
for i=1:2:nargs
    switch args{i},
        case 'pa_limit',   pa_limit=args{i+1};
        case 'hprior', Hprior = args{i+1};
        case 'sprior', Sprior = args{i+1};
    end
end


if isempty(Hprior)
    Hprior = ones(nGenes);
    Hprior = setdiag(Hprior,0);
end
%% precompute local score


LS=cell(size(expdata,1),1);

dummy=uint16(nGenes+1); % set dummy number for no-parent
nopavec=dummy:dummy+pa_limit-1;
padim=double(repmat(dummy,1,pa_limit));

X1=ones(nSamples,1);
for i=1:size(expdata,1)
    % i
    parents_cand = find(Hprior(:,i));
    LS{i}.score = []; % local score
    LS{i}.parent = []; % parent combination

    Spvec=Sprior(:,i)';
    y = expdata(i,:)';
    for j=0:min([pa_limit,length(parents_cand)])
        
        pasets_cand = VChooseK(parents_cand,j); % get parent combination
        
        % set prior into localscore
        if j==1
            localscore = Spvec(pasets_cand)';
        else
            localscore = sum(Spvec(pasets_cand),2);
        end
        
        % calc local score
        if isempty(pasets_cand)
            
            
            beta = X1 \ y;
            
            yhat = X1 *beta;
            residuals = y - yhat;
            params(2)=sqrt(residuals'*residuals/nSamples);
            
            loglik = -normlike(params,residuals);
            localscore = loglik - 0.5*length(beta)*log(nSamples);
        else
            for k=1:size(pasets_cand,1)

                X = [ones(1,nSamples);expdata(pasets_cand(k,:),:)]';
                beta = X \ y;

                yhat = X*beta;
                residuals = y - yhat;                
                loglik = -my_normlike(residuals,0,sqrt(residuals'*residuals/nSamples));
                localscore(k) = localscore(k) + loglik - 0.5*length(beta)*log(nSamples);% should be replaced by my_normlike
            end
        end
        
        LS{i}.score = [LS{i}.score;localscore];
        patemp=[pasets_cand,repmat(nopavec(end-pa_limit+j+1:end),max([size(pasets_cand,1),1]),1)];
        LS{i}.parent = [LS{i}.parent;patemp];

    end
    max_LL=max(LS{i}.score);
    LS{i}.score=LS{i}.score-max_LL;
    LS{i}.max_LL = max_LL;
    
end


% make lookup table and index and sort
for i=1:length(LS)
    locmat = false(length(LS{i}.score),length(LS)+pa_limit);
    for j=1:length(LS)+pa_limit
        locmat(:,j)=any(LS{i}.parent==j,2);
    end
    locmat = sparse(locmat);
    LS{i}.lookup = locmat;
    LS{i}.parent(LS{i}.parent>nopavec(1))=nopavec(1);% re-pudding blanks with the same value

    % make index
    temp = cell(pa_limit,1);
    for k = 1:pa_limit
        temp{k}=LS{i}.parent(:,k);
    end
    LS{i}.indx = double(sub2ind(padim,temp{:}));
    
    
    %sort
    [~, bi]=sort(LS{i}.indx);
    LS{i}.indx=LS{i}.indx(bi);
    LS{i}.parent=LS{i}.parent(bi,:);
    LS{i}.score= LS{i}.score(bi);
    LS{i}.lookup = LS{i}.lookup(bi,:);
end

end


