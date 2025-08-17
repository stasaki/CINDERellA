function [r]=getpermTh(expdata,snps,genepos,snpspos,cisDist,pvOutputThreshold)
%getpermTh Summary of this function goes here
% INPUTS:
%	expdata:	parameter object
%   snps: save result?
%   genepos: save result?
%   snpspos: save result?
%   cisDist: save result?
%   pvOutputThreshold: save result?
% OUTPUTS:
%	snpind:	r
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License

% load expdata and snpdata into SlicedData object

Gene = SlicedData;
Gene.fileDelimiter = char(9); % the TAB character
Gene.fileOmitCharacters = NaN; % denote missing values;
Gene.fileSkipRows = 0; % one row of column labels
Gene.fileSkipColumns = 0; % one column of row labels
Gene.fileSliceSize = 10000; % read file in pieces of 2,000 rows
Gene.LoadBFile( expdata,genepos.symbol );


Snp = SlicedData;
Snp.fileDelimiter = char(9); % the TAB character
Snp.fileOmitCharacters = NaN; % denote missing values;
Snp.fileSkipRows = 0; % one row of column labels
Snp.fileSkipColumns = 0; % one column of row labels
Snp.fileSliceSize = 10000; % read file in pieces of 2,000 rows
Snp.LoadBFile( snps,snpspos.name );

cvrt = SlicedData;
%

nPerm=1000;

tic;

geneperm = Gene.Clone();
permR=Matrix_eQTL_engine_cis_maxR_perm( Snp, ...
    geneperm, ...
    cvrt, ...
    snpspos, ...
    genepos, ...
    cisDist,...
    nPerm);

toc;
if size(permR,1)~=1
    permR=nanmax(permR);
end

permR=sort(permR);
r=permR(end-round(pvOutputThreshold*nPerm));


