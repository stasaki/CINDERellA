function DcalcLSC_cmat( Param )
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

load(Param.savepath.simdata);

num_ls = 0;
for i=1:size(Hprior.skeleton,1)
    num_ls = num_ls+ my_combin(sum(Hprior.skeleton(:,i)),Param.pa_limit);
end
% max_ls = my_combin(Param.nGene-1,Param.pa_limit)*Param.nGene;




if (isdeployed)
    cd('/scratch/stasaki/');
    % cd('/scratch/sgeadmin/');
    workdir =  strrep(Param.savepath.ls,'.mat','');
    mkdir(workdir);
    cd(workdir);
    
    disp('Exporting Exppression data as a text file...')
    dlmwrite('./exp.txt', expdata.data, 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    fprintf('... Done!\n\n');

    disp('Exporting number of state as a text file...')
    dlmwrite('./ns.txt', expdata.ns, 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    fprintf('... Done!\n\n');
    
    disp('Exporting Hard prior as a text file...')
    dlmwrite('./hp.txt', full(Hprior.skeleton), 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    clear Hprior
    fprintf('... Done!\n\n');
    
    status = unix(['DcalcLS ',num2str(Param.nGene),' ',...
        num2str(Param.nSample),' exp.txt ',...
        'hp.txt ',spfile_loc,' ','ns.txt ',...
        num2str(Param.pa_limit),' 1 ',num2str(num_ls)]);
else
    workdir = strrep(Param.savepath.ls,'.mat','');
    mkdir(workdir);
    cd(workdir);
    disp('Exporting Exppression data as a text file...')
    dlmwrite('./exp.txt', expdata.data, 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    fprintf('... Done!\n\n');
    
    disp('Exporting number of state as a text file...')
    dlmwrite('./ns.txt', expdata.ns, 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    fprintf('... Done!\n\n');
    
    disp('Exporting Hard prior as a text file...')
    dlmwrite('./hp.txt', full(Hprior.skeleton), 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    clear Hprior
    fprintf('... Done!\n\n');
    
    if  ~strcmp('[]',spfile_loc)
    disp('Exporting Soft prior as a text file...')
    dlmwrite('./sp.txt', full(Prior), 'delimiter', '\t', ...
        'precision', 6,'newline','unix');
    clear Prior
    end
    
    fprintf('... Done!\n\n');    
    status = unix(['/home/shinya/Resource/compbiotools/expression/bayesian_network/cbin/linux/DcalcLS ',num2str(Param.nGene),' ',...
        num2str(Param.nSample),' exp.txt ',...
        'hp.txt ',spfile_loc,' ','ns.txt ',...
        num2str(Param.pa_limit),' 1 ',num2str(num_ls)]);
end




if status~=0
    disp('error in eSNPnet')
    exit;
end



delete('exp.txt')
delete('hp.txt')

cd(currentFolder);


end

