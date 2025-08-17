function CcalcLSC_cmat( Param,expdata )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

currentFolder = pwd;
if strcmp('no',Param.hprior)
    Hprior.skeleton = ones(Param.nGene);
    Hprior.skeleton=setdiag(Hprior.skeleton,0);
else
    load(Param.savepath.hprior);
end

spfile_loc='xxx';


%load(Param.savepath.simdata);

num_ls = 0;
for i=1:size(Hprior.skeleton,1)
    num_ls = num_ls+ my_combin(sum(Hprior.skeleton(:,i)),Param.pa_limit);
end
% max_ls = my_combin(Param.nGene-1,Param.pa_limit)*Param.nGene;



% 
% if (isdeployed)
%     cd('/scratch/stasaki/');
%     % cd('/scratch/sgeadmin/');
%     workdir = [Param.id,'_LS_',Param.netid,'_hp',Param.hprior, num2str(randi(10000000,1))];
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
%     status = unix(['ClinP ',num2str(Param.nGene),' ',...
%         num2str(Param.nSample),' exp.txt ',...
%         'hp.txt ',spfile_loc,' ',...
%         num2str(Param.pa_limit),' 1 ',num2str(num_ls)]);
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
    % should set path to binary 
    %/Resource/compbiotools/expression/bayesian_network/cbin/linux
    status = unix(['ClinP ',num2str(Param.nGene),' ',...
        num2str(Param.nSample),' exp.txt ',...
        'hp.txt ',spfile_loc,' ',...
        num2str(Param.pa_limit),' 1 ',num2str(num_ls)]);
% end


delete('exp.txt')
delete('hp.txt')

if status~=0
    disp('error')
else
    unix('touch done.txt')
end
cd(currentFolder);

end

