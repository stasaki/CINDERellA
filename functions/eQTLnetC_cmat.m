function eQTLnetC_cmat( Param )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

currentFolder = pwd;
if strcmp('no',Param.hprior)
    Hprior.skeleton = ones(Param.nGene);
    Hprior.skeleton=setdiag(Hprior.skeleton,0);
else
    load(Param.savepath.hprior);
end

if isfield(Param.savepath,'sprior')
    load( Param.savepath.sprior);
    spfile_loc = 'sp.txt';
else
    spfile_loc='[]';
end

load(Param.savepath.esnp)
load(Param.savepath.simdata)


num_ls = 0;
for i=1:size(Hprior.skeleton,1)
    num_ls = num_ls+ my_combin(sum(Hprior.skeleton(:,i)),Param.pa_limit);
end

nSnp = num2str(length(find(sum(SNP.snpind)>0)));
% 
% if (isdeployed)
%     cd('/scratch/stasaki/');
%     % cd('/scratch/sgeadmin/');
%     workdir = [date_string(),'LS_',Param.netid,'_hp',num2str(Param.hprior), num2str(randi(100000000,1))];
%     mkdir(workdir);
%     cd(workdir);
%     
%     disp('Exporting Exppression data as a text file...')
%     dlmwrite('./exp.txt', expdata.data, 'delimiter', '\t', ...
%         'precision', 6,'newline','unix');
%     fprintf('... Done!\n\n');
%     
%     disp('Exporting Hard prior as a text file...')
%     dlmwrite('./hp.txt', full(Hprior.skeleton), 'delimiter', '\t', ...
%         'precision', 6,'newline','unix');
%     clear Hprior
%     fprintf('... Done!\n\n');
%     
%     disp('exporting Genotype data as a text file...')    
%     ext_indx= sum(SNP.snpind)>0;
%     dlmwrite('./snp.txt', snps.data(ext_indx,:), 'delimiter', '\t', ...
%         'precision', 6,'newline','unix');
%     [geneid snpid]=find(SNP.snpind(:,ext_indx));
%     dlmwrite('./snp2gene.txt', [geneid snpid], 'delimiter', '\t', ...
%         'precision', 6,'newline','unix');
%     fprintf('... Done!\n\n');
% 
%     status = unix(['eSNPnet ',num2str(Param.nGene),' ',nSnp,' ',...
%         num2str(Param.nSample),' exp.txt snp.txt ',...
%         'snp2gene.txt hp.txt ',spfile_loc,' ',...
%         num2str(Param.pa_limit),' ',num2str(SNP.lrThresh(1))]);
% else
    workdir = strrep(Param.savepath.ls,'.mat','');
    mkdir(workdir);
    cd(workdir);
    
    disp('Exporting Exppression data as a text file...')
    dlmwrite('./exp.txt', expdata.data, 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    fprintf('... Done!\n\n');
    
    disp('Exporting Hard prior as a text file...')
    dlmwrite('./hp.txt', full(Hprior.skeleton), 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    clear Hprior
    fprintf('... Done!\n\n');
    
    
    disp('exporting Genotype data as a text file...')
    ext_indx= sum(SNP.snpind)>0;
    dlmwrite('./snp.txt', snps.data(ext_indx,:), 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    [geneid, snpid]=find(SNP.snpind(:,ext_indx));
    dlmwrite('./snp2gene.txt', [geneid snpid], 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    fprintf('... Done!\n\n');
    
    if  ~strcmp('[]',spfile_loc)
        disp('Exporting Soft prior as a text file...')
        dlmwrite('./sp.txt', full(Prior), 'delimiter', '\t', ...
            'precision', 6,'newline','unix');
        clear Prior
    end
    % should set path to binary 
    %/Resource/compbiotools/expression/bayesian_network/cbin/linux
    
    status = unix(['eSNPnetP ',num2str(Param.nGene),' ',nSnp,' ',...
        num2str(Param.nSample),' exp.txt snp.txt ',...
        'snp2gene.txt hp.txt ',spfile_loc,' ',...
        num2str(Param.pa_limit),' ',num2str(SNP.lrThresh(1)),' 1 ',num2str(num_ls)]);
% end



if status~=0
    disp('error in eSNPnet')
    exit;
end

delete('exp.txt')
delete('hp.txt')
delete('snp.txt')
unix('touch done.txt')
cd(currentFolder);



end

