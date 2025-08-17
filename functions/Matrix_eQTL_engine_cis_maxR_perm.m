function [nullRs]=Matrix_eQTL_engine_cis_maxR_perm(snps,gene,cvrt,snpspos,genepos,cisDist,nPerm)
% INPUTS:
%	snps:	SlicedData object for genotype matrix
%	gene: 	SlicedData object for expression matrix
%   cvrt:   SlicedData object for covariate matrix
%   snpspos:    chromosome position of snps
%   genepos:    chromosome position of genes
%   cisDist:    cis-SNPs range
%   nPerm:      the number of permutation
% OUTPUTS:
%	nullRs:		null distribution of SNP-gene correlation
%
% The Original code from Matrix eQTL by Andrey A. Shabalin
% http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
% Build 1.2.0
%
% The original script has been modified to get null R distribution by
% permuting the samples.
% Modification made on Nov. 2012 by:
%
% Shinya Tasaki, Ph.D.



% Set missing parameters

verbose = false;
errorCovariance = [];


gene = gene.Clone();
%snps = snps.Clone();
cvrt = cvrt.Clone();


snps_process = @impute_row_mean;

if(verbose)
	tic
	status = @status_T;
else
	status = @status_F;
end;

% Check dimensions
status('Checking input data dimensions',[]);
if(snps.nCols()*snps.nRows() == 0)
	error('Empty genotype dataset');
end;
if(gene.nCols()*gene.nRows() == 0)
	error('Empty expression dataset');
end;
if(snps.nCols ~= gene.nCols)
	error('Different number of samples in the genotype and gene expression files');
end;
if(cvrt.nRows>0)
	if(snps.nCols ~= cvrt.nCols)
		error('Wrong number of samples in the file with covariates');
	end;	
end;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% error covariance processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~isempty(errorCovariance))
	status('Processing the errorCovariance matrix');
	if(size(errorCovariance,1)~=size(errorCovariance,2))
		error('The covariance matrix is not square');
	end;
	if(size(errorCovariance,1)~=snps.nCols)
		error('The covariance matrix size does not match the data');
	end;	
	% test for symmetry
	if(~all(all(errorCovariance==errorCovariance')))
		error('The covariance matrix is not symmetric');
	end;
	[v,d] = eig(errorCovariance);	
	%  errorCovariance == v*d*v'
	%  errorCovariance^0.5 == v*sqrt(d)*v'
	%  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v'
	d = diag(d);
	if(any(d<=0))
		error('The covariance matrix is not positive definite');
	end;
	correctionMatrix = v*diag(1./sqrt(d))*v';
	clear v d;
else 
	clear correctionMatrix;
	correctionMatrix = [];
end;

%%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% covariates processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add constant as a covariate
cvrt.SetNanRowMean();
cvrt.CombineInOneSlice(); %
cvrt = [ones(1,snps.nCols);cvrt.dataSlices{:}];

% Correct for the error covariance structure
if(~isempty(correctionMatrix))
	status('Rotating cvrt based on the errorCovariance matrix');
	cvrt = cvrt * correctionMatrix;
end;

% Orthonormalize covariates
status('Orthonormalizing covariates');
[q, r] = qr(cvrt',0);
if(min(abs(diag(r))) < eps(class(r))*snps.nCols())
	error('Colinear or zero covariates detected.');
end;
cvrt = q';
clear q;

%%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% gene expression processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Impute gene expression
status('Imputing missing expression');
gene.SetNanRowMean();
% Correct for the error covariance structure
if(~isempty(correctionMatrix))
	status('Rotating expression based on the errorCovariance matrix');
	gene.RowMatrixMultiply(correctionMatrix);
end;

gene.RowStandardizeCentered();

% Orthogonolize expression w.r.t. the covariates
status('Orthogonolizing expression w.r.t. covariates');
for sl = 1:gene.nSlices
	slice = gene.dataSlices{sl};
	slice = slice - (slice*cvrt')*cvrt;
	gene.dataSlices{sl} = slice;
end;
clear sl d slice;
status('Standardizing expression');
gene.RowRemoveZeroEps();
gene.RowStandardizeCentered();


%%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matching gene and SNPs locations    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

status('Matching data files and location files')
	
% names in the input data
gene_names = vertcat(gene.rowNameSlices{:});
snps_names = vertcat(snps.rowNameSlices{:});

% match with the location data
usedgene = false(size(genepos.symbol,1),1);
[ai, genematch2, genematch] = intersect(gene_names, genepos.symbol);
usedgene( genematch ) = 1;
if isempty(genematch)
    error('Gene names do not match those in the gene location file.');
end
% disp([num2str(length(genematch)),' of ',num2str(length(gene_names)),' genes matched']);
	
usedsnps = false(size(snpspos.name,1),1);
[ai, snpsmatch2, snpsmatch] = intersect(snps_names, snpspos.name);
usedsnps( snpsmatch ) = 1;
if isempty(snpsmatch)
    error('SNP names do not match those in the SNP location file.');
end
% disp([num2str(length(snpsmatch)),' of ',num2str(length(snps_names)),' snps matched']);
	
% find accessed chr names
chrNames = unique([unique(snpspos.chrm_snp(usedsnps));unique(genepos.chrm_probe(usedgene))]);

% match chr names
genechr=zeros(length(genepos.chrm_probe),1);
snpschr=zeros(length(snpspos.chrm_snp),1);
for i=1:length(chrNames)
    [ai, bi, ci] = my_intersect(genepos.chrm_probe,chrNames(i));
    genechr(bi)=i;
    [ai, bi, ci] = my_intersect(snpspos.chrm_snp,chrNames(i));
    snpschr(bi)=i;
end
	


% max length of chromosome
chrMax = max( [snpspos.pos(usedsnps); genepos.txEnd(usedgene)]);
	
% set single number location
genepos2 = [genepos.txStart,genepos.txEnd];
genepos2 = genepos2 + repmat((genechr-1)*(chrMax+max(cisDist)),1,2);

snpspos2 = snpspos.pos;
snpspos2 = snpspos2 + (snpschr-1)*(chrMax+max(cisDist));

snps_pos = zeros(length(snps_names),1);
snps_pos(snpsmatch2) = snpspos2(snpsmatch);
snps_pos(snps_pos==0) = (length(chrNames)+1)*(chrMax+max(cisDist));
	
gene_pos = zeros(length(gene_names),2);
gene_pos(genematch2,:) = genepos2(genematch,:);
gene_pos(gene_pos(:,1)==0,:) = (length(chrNames)+2)*(chrMax+max(cisDist));


gene_strand = ones(length(gene_names),1);
gene_strand(genematch2,:) = genepos.strand(genematch,:);

gene_pos(gene_strand==-1,:)=gene_pos(gene_strand==-1,[2,1]);

clear genematch genematch2 usedgene snpsmatch snpsmatch2 usedsnps chrNames genechr snpschr chrMax genepos2 snpspos2
	
% Slice it back.
geneloc = cell(gene.nSlices,1);
genestrand = cell(gene.nSlices,1);
gene_offset = 0;
for gc=1:gene.nSlices
    nr = length(gene.rowNameSlices{gc});
    geneloc{gc} = gene_pos(gene_offset + (1:nr),:);
    genestrand{gc} = gene_strand(gene_offset + (1:nr),:);
    gene_offset = gene_offset + nr;
end
clear gc gene_offset gene_pos

snpsloc = cell(snps.nSlices,1);
snps_offset = 0;
for sc = 1:snps.nSlices
    nr = length(snps.rowNameSlices{sc});
    snpsloc{sc} = snps_pos(snps_offset + (1:nr), :);
    snps_offset = snps_offset + nr;
end
clear sc snps_offset snps_pos

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

status('Performing eQTL analysis.');
snps_offset = 0;
nSample = snps.nCols;
nullRs = nan(snps.nSlices*gene.nSlices,nPerm);
idx=1;
for sc = 1:snps.nSlices
    gene_offset = 0;
    
    % prepare snps / dummies
    cursnps = snps_process( snps.dataSlices{sc} );
    for d = 1:length(cursnps)
        if(~isempty(correctionMatrix))
            cursnps{d} = cursnps{d} * correctionMatrix;
        end;
        cursnps{d} = cursnps{d} - (cursnps{d}*cvrt')*cvrt;
        for w = 1:(d-1)
            cursnps{d} = cursnps{d} - ...
                bsxfun(@times, sum(cursnps{d}.*cursnps{w},2), cursnps{w});
        end;
        cursnps{d} = RowStandardizeCentered(cursnps{d});
    end;
    
    nrcs = size(cursnps{1},1);
    
    for gc = 1:gene.nSlices
        curgene = gene.dataSlices{gc};
        nrcg = size(curgene,1);
        
        srep = repmat( snpsloc{sc},1, nrcg);
        iscis = ( bsxfun(@times,srep,genestrand{gc}(:,1)') > bsxfun(@times,repmat((geneloc{gc}(:,1) - bsxfun(@times,genestrand{gc}(:,1),cisDist(1)))',nrcs,1),genestrand{gc}(:,1)')) &...
            ( bsxfun(@times,srep,genestrand{gc}(:,1)') < bsxfun(@times,repmat((geneloc{gc}(:,2) + bsxfun(@times,genestrand{gc}(:,1),cisDist(2)))',nrcs,1),genestrand{gc}(:,1)'));
        clear srep
        
        if isempty(find(iscis, 1))
            
        else

                for np=1:nPerm
                    %np
                cor = (cursnps{1}*curgene(:,randperm(nSample))');
                nullRs(idx,np)=max(abs(cor(iscis)));
                
                end

        end
        idx=idx+1;
        
        gene_offset = gene_offset + nrcg;
    end;
    snps_offset = snps_offset + nrcs;
end;

clear cor r2 select_cis select_tra test pv cursnps curgene

clear cor select adump pdump pv testmatrix

end

function status_T(text,~)
	if(nargin==1)
		disp(['Task	finished in ' num2str(toc) ' seconds']);
	end;
	disp(text);
	tic;
end

function status_F(~,~)
end

function y = impute_row_mean(x)
	if(any(x))
		rowmean = nanmean(x,2);
		rowmean(isnan(rowmean)) = 0;
		for j=1:size(x,2)
			where1 = isnan(x(:,j));
			x(where1,j) = rowmean(where1);
		end
	end;
	y = {x};
end

function y = RowStandardizeCentered(x)
	div = sqrt(sum(x.^2,2));
	div(div==0) = 1;
	y = bsxfun(@rdivide, x, div);
end

