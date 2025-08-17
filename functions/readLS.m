function [ LS ] = readLS( Param )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

currentFolder = pwd;

workdir =  strrep(Param.savepath.ls,'.mat','');

cd(workdir);

nGenes=Param.nGene;
pa_limit = Param.pa_limit;
dummy=uint16(nGenes+1); % set dummy number for no-parent
nopavec=dummy:dummy+pa_limit-1;
padim=double(repmat(dummy,1,pa_limit));

LS=cell(nGenes,1);

for i=1:nGenes
    
    fid = fopen(['LS',num2str(i),'.cmat']);
    LS{i}.score = fread(fid,inf,'double');
    fclose(fid);
    
    
    fid = fopen(['Pa',num2str(i),'.cmat']);
    LS{i}.parent = fread(fid,[length(LS{i}.score),pa_limit],'int');
    fclose(fid);
    
    
    max_LL=max(LS{i}.score);
    LS{i}.score=LS{i}.score-max_LL;
    LS{i}.max_LL = max_LL;
    
end

% make lookup table
for i=1:length(LS)
    locmat = false(length(LS{i}.score),length(LS)+pa_limit);
    for j=1:length(LS)+pa_limit
        locmat(:,j)=any(LS{i}.parent==j,2);
    end
    locmat = sparse(locmat);
    LS{i}.lookup = locmat;
    LS{i}.parent(LS{i}.parent>nopavec(1))=nopavec(1);
    %LS{i}.parent=uint16(LS{i}.parent);
    
    % make index
    temp = cell(pa_limit,1);
    for k = 1:pa_limit
        temp{k}=LS{i}.parent(:,k);
    end
    LS{i}.indx = double(sub2ind(padim,temp{:}));
    
    
    
    % make index
    %LS{i}.indx = subv2indMinka(padim,double(LS{i}.parent));
    %sort
    [~, bi]=sort(LS{i}.indx);
    LS{i}.indx=LS{i}.indx(bi);
    LS{i}.parent=LS{i}.parent(bi,:);
    LS{i}.score= LS{i}.score(bi);
    LS{i}.lookup = LS{i}.lookup(bi,:);
end

cd(currentFolder);

end

