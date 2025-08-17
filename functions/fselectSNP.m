function [ snpind,rThresh,lrThresh,corvec ] = fselectSNP( expdata,genodata,varargin)
%fselectSNP Summary of this function goes here
% INPUTS:
%	expdata:	parameter object
%   genodata: save result?
%   pvOutputThreshold: save result?
%   genepos: save result?
%   snpspos: save result?
%   cisDist: save result?
%   method: save result?
%   perm: save result?
% OUTPUTS:
%	snpind:	mcmc result
%	rThresh:	mcmc result
%	lrThresh:	mcmc result
%	corvec:	mcmc result
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License


% rng('shuffle')
[nGenes, nSamples] = size(expdata);

method='all';
pa_limit=3;
pvOutputThreshold = 0.01;
cis=1;
perm = 0;

args = varargin;
nargs = length(args);


for i=1:2:nargs
    switch args{i},
        case 'pvOutputThreshold',   pvOutputThreshold=args{i+1};
        case 'genepos', genepos = args{i+1};
        case 'snpspos', snpspos = args{i+1};
        case 'cisDist', cisDist=args{i+1};
        case 'method', method = args{i+1};
        case 'perm' ,perm=args{i+1};
    end
end



%% precompute local score

nVarTested =1;
rThresh = zeros(pa_limit+1,1);
lrThresh = zeros(pa_limit+1,1);


if perm ==1
    [r]=getpermTh(expdata,genodata,genepos,snpspos,cisDist,pvOutputThreshold);
    rThresh(:)=r;
    lrThresh(:)=-nSamples*log(1-r.^2)/2;
    clear r;
else
    for j=0:pa_limit
        nCov = j;
        dfFull = nSamples - nCov - nVarTested;
        tThresh = -tinv(pvOutputThreshold/2,dfFull);
        rThresh(j+1) = sqrt(  tThresh.^2 ./  (dfFull + tThresh.^2)  );
        clear tThresh;
    end
    lrThresh(:)=-nSamples*log(1-rThresh.^2)/2;
end

if strcmp('pca',method)
   pcs = cell(nGenes,1);
   numpcs = zeros(nGenes,1);
else
    snpind=sparse(false(nGenes,size(genodata,1)));
end

genotest=genodata;
corvec = nan(size(expdata,1),1);
for i=1:size(expdata,1)
    %i
    genetest = expdata(i,:);

    
    
    if cis ==1
        genepos2 = [genepos.txStart(i),genepos.txEnd(i)];
        if genepos.strand(i)==-1
            genepos2 = genepos2([2,1]);
        end
        genepos2 = genepos2+cisDist.*[-1,1].*double(genepos.strand(i));
        snpsflag = snpspos.chrm_snp==genepos.chrm_probe(i)&snpspos.pos>min(genepos2)&snpspos.pos<max(genepos2);
        genotest = genodata(snpsflag,:);
    
        
        if ~isempty(genotest)
            switch method
                case 'all'
                    
                    snpind(i,snpsflag)=1;
                    
                case 'pca'
                    if size(genotest,1)~=1 %% need to modify when only one snps
                        [~,SCORE,latent] = princomp(zscore(genotest'));
                        genotest=SCORE(:,cumsum(latent)./sum(latent) < 0.9)';
                    end
                    pcs{i}=genotest;
                    numpcs(i)=size(genotest,1);
                case 'LD'
                    genosnppos = snpspos.pos(snpsflag);
                    [genosnppos bi]=sort(genosnppos);
                    genotest=genotest(bi,:);

                    [r p]=corr(genotest');
                    
                    ldtemp = zeros(length(r));
                    for k=1:length(r)
                        ldtemp(k,:)=p(k,:)>0.01;
                        ldtemp(k,k)=0;
                        temp=find(ldtemp(k,1:k-1),1,'last');
                        ldtemp(k,1:temp)=1;
                        temp=find(ldtemp(k,k+1:end),1,'first');
                        ldtemp(k,k+1+temp:end)=1;
                        
                    end
                    
                    
                    ldtemp=~ldtemp;
                    ld = ones(1,length(r));
                    s=1;
                    for k=1:length(r)-1
                        if any((ldtemp(:,k)+ldtemp(:,k+1))==2)
                            ld(k+1)=s;
                        else
                            s=s+1;
                        end
                    end
                    
                    
                    maf=sum(genotest,2)./(size(genotest,2)*2);
                    maf(maf>0.5)=1-maf(maf>0.5);
                    genoind=false(length(r),1);
                    for k=unique(ld)
                        maxmaf=max(maf(ld == k));
                        indtemp=find((ld == k)'&(maf==maxmaf));
                        if ~isempty(indtemp)
                            genoind(indtemp(1))=1;
                        end
                    end
                    snpsind=find(snpsflag);
                    snpind(i,snpsind(genoind))=1;
                    
                    
                case 'no-parent-max'
                    [q, r] = qr(ones(size(expdata,2),1),0);
                    
                    
                    
                    genetemp = genetest - (genetest*q)*q';
                    genetemp = RowStandardizeCentered(genetemp);
                    snptest = genotest - (genotest*q)*q';
                    snptest = RowStandardizeCentered(snptest);
                    
                    
                    cor = (snptest*genetemp');
                    
                    [ai, bi]=max(abs(cor));

                    snpsind=find(snpsflag);
                    snpind(i,snpsind(bi))=1;
                    corvec(i)=ai;
                   
 
            end
        end
    end
    
end


if strcmp('pca',method)
    
end

end

