function [LS, LS_null] = eQTLnet(expdata, genodata, snpind,rThresh,varargin)
% eQTLnet pre-computing local score for network fragments with genotype
% INPUTS:
%	expdata:	gene expression data
%   genodata: genotype matrix
%   snpind: snp to gene relations
%   rThresh: R threshold
%   pa_limit: the number of parents allowed
%	hprior:     matrix of acceptable edges: 1-allowed, 0-banned
%   sprior:   log of prior for the egdge
% OUTPUTS:
%	LS:		local score object with snp
%   LS_null: local score object without snp
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License

% rng('shuffle')
[nGenes, nSamples] = size(expdata);


pa_limit=3;
Sprior = zeros(nGenes);
Hprior = [];

args = varargin;
nargs = length(args);


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


%% common settings
LS=cell(size(expdata,1),1);
LS_null=LS;
dummy=uint16(nGenes+1); % set dummy number for no-parent
nopavec=dummy:dummy+pa_limit-1;
padim=double(repmat(dummy,1,pa_limit));

%% precompute local score

cvrtmp1 = ones(nSamples,1);
[q1, r]= qr(cvrtmp1,0);
for i=1:size(expdata,1)
    i
    parents_cand = find(Hprior(:,i));
    
    LS{i}.score = [];
    LS{i}.parent = [];
    %     LS{i}.indx = [];   % linear index of parent combination
    LS{i}.snps = [];
    genetest = expdata(i,:)';
    
    snploc=find(snpind(i,:));
    Spvec=Sprior(:,i)';
    for j=0:min([pa_limit,length(parents_cand)])
        
        pasets_cand = VChooseK(parents_cand,j);
        
        if j==1
            localscore = Spvec(pasets_cand)';
        else
            localscore = sum(Spvec(pasets_cand),2);
        end
        
        selectedpcs = zeros(max([1,size(pasets_cand,1)]),1);
        if isempty(pasets_cand)
            
            if ~isempty(snploc)
                
                % Orthonormalize covariates
                genetemp = genetest' - (genetest'*q1)*q1';
                genetemp = RowStandardizeCentered(genetemp);
                snptest = genodata(snploc,:) - (genodata(snploc,:)*q1)*q1';
                snptest = RowStandardizeCentered(snptest);
                cor = (snptest*genetemp');
                
                [ai, bi]=max(abs(cor));
                if ai >=  rThresh
                    selectedpcs(1) = snploc(bi);
                else
                    bi=[];
                end
                
                X = [cvrtmp1,genodata(snploc(bi),:)'];
                beta = X \ genetest;
                yhat = X*beta;
                residuals = genetest - yhat;
                loglik = -my_normlike(residuals,0,sqrt(residuals'*residuals/nSamples));
                localscore =  loglik - 0.5*length(beta)*log(nSamples);
            else
                beta = cvrtmp1 \ genetest;
                yhat = cvrtmp1*beta;
                residuals = genetest - yhat;
                loglik = -my_normlike(residuals,0,sqrt(residuals'*residuals/nSamples));
                localscore = loglik - 0.5*length(beta)*log(nSamples);
            end
        else
            for k=1:size(pasets_cand,1)
                cvrtmp = [cvrtmp1,expdata(pasets_cand(k,:),:)'];
                if ~isempty(snploc)
                    
                    % Orthonormalize covariates
                    [q, r]= qr(cvrtmp,0);
                    genetemp = genetest' - (genetest'*q)*q';
                    genetemp = RowStandardizeCentered(genetemp);
                    snptest = genodata(snploc,:) - (genodata(snploc,:)*q)*q';
                    snptest = RowStandardizeCentered(snptest);
                    cor = (snptest*genetemp');
                    [ai, bi]=max(abs(cor));
                    if ai >=  rThresh
                        selectedpcs(k) = snploc(bi);
                    else
                        bi=[];
                    end
                    
                    
                    X = [cvrtmp,genodata(snploc(bi),:)'];
                    beta = X \ genetest;
                    yhat = X*beta;
                    residuals = genetest - yhat;
                    loglik = -my_normlike(residuals,0,sqrt(residuals'*residuals/nSamples));
                    localscore(k) = localscore(k)+ loglik - 0.5*length(beta)*log(nSamples);
                    
                else
                    beta = cvrtmp \ genetest;
                    yhat = cvrtmp*beta;
                    residuals = genetest - yhat;
                    loglik = -my_normlike(residuals,0,sqrt(residuals'*residuals/nSamples));
                    localscore(k) = localscore(k)+ loglik - 0.5*length(beta)*log(nSamples);
                end
            end
        end
        
        
        LS{i}.score = [LS{i}.score;localscore];
        patemp=[pasets_cand,repmat(nopavec(end-pa_limit+j+1:end),max([size(pasets_cand,1),1]),1)];

        LS{i}.parent = [LS{i}.parent;patemp];
        LS{i}.snps = [LS{i}.snps;selectedpcs];
    end
    
    % normalize
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
    %     LS{i}.parent = single(LS{i}.parent);
    
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
    LS{i}.snps = LS{i}.snps(bi,:);
end



end



