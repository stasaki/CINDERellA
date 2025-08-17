function [sampled_graphs, LL_retain] = learn_struct_mcmc_gibbs_block4_con(varargin)
%learn_struct_mcmc_gibbs_block4_con learning network structure by mcmc
% INPUTS:
%	nsamples:	the number of samples collected from the chain
%   init_dag: initail network structure
%   LS: local score object
%   hprior: matrix of acceptable edges: 1-allowed, 0-banned
%   pa_limit: the number of maximum parents
%   maxlen: the maximum number of linear connected nodes to be modified (3-c4PB,2-c3PB)
%   runtime: running time
%   sample_duration: time intervel collecting samples
% OUTPUTS
%	sampled_graphs:	sampled networks
%	LL_retain:	log likelihood of sampled networks
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License



% debug=0;

args = varargin;
nargs = length(args);

% set default parameters
Hprior=[];
runt=inf;

% set paramters
for i=1:2:nargs
    switch args{i},
        case 'nsamples',   nsamples = args{i+1}; %nsamples=100
        case 'init_dag',   dag = args{i+1};
        case 'LS',       LS = args{i+1};
        case 'hprior', Hprior = args{i+1};
        case 'pa_limit',    pa_limit = args{i+1};
        case 'maxlen', maxlen=args{i+1};
        case 'runtime', runt=args{i+1};
        case 'sample_duration', sample_duration=args{i+1};
    end
end




% set number of iteration
sampled_graphs = cell(1, nsamples);
dag = double(dag);



% calc reachability of inithial network- counting number of path from nodeX
% to nodeY
n=length(dag); % get number of nodes in network
[ai, bi] = find(dag==1);
cmat = zeros(n); % cmat stores number of path between nodes
cmat(diag(ones(n,1))==1)=1;
for i=1:length(ai)
    cmat=cmat+cmat(:,ai(i))*cmat(bi(i),:);
end


% hard prior
if isempty(Hprior)
    Hprior = ones(n);
    Hprior = setdiag(Hprior,0);
end


LL_retain = -inf(nsamples,1);
dummy=n+1; % specify dummy parent indx for no-parent
padim=repmat(dummy,1,pa_limit);


% calc score for initial network
LL = score_dag(dag,LS,pa_limit,dummy,padim);


% compute sample space code
load('net4code.mat');
load('net3code.mat');

% pre declaration
n1_3indx = cell(4,1);
n2_3indx = cell(4,1);
n3_3indx = cell(4,1);
n1_3LS_indx = cell(4,1);
n2_3LS_indx = cell(4,1);
n3_3LS_indx = cell(4,1);
pb_3n1 = cell(25,1);
pb_3n2 = cell(12,1);

n1_4indx = cell(8,1);
n2_4indx = cell(8,1);
n3_4indx = cell(8,1);
n4_4indx = cell(8,1);
n1_4LS_indx = cell(8,1);
n2_4LS_indx = cell(8,1);
n3_4LS_indx = cell(8,1);
n4_4LS_indx = cell(8,1);
pb_4n1 = cell(543,1);
pb_4n2 = cell(200,1);
pb_4n3 = cell(48,1);

temp2=repmat(dummy,1,pa_limit);
sampind=1;


node_temp=nan(maxlen,1);
next_sample_time = sampind*sample_duration;
tic
while 1
    
    
    % all find path
    [ai, bi]=find(cmat>0);
    % pick one path
    indx = randi(length(ai),1,1);
    from=ai(indx);
    to=bi(indx);
    
    k=1;
    node_temp(:)=nan;
    haspath = cmat(:,to)'>0;
    while from ~= to
        if  k >= maxlen+1
            break;
        end
        node_temp(k)=from;
        next_nodes = find(dag(from,:)&haspath);
        len_next_nodes = length(next_nodes);
        if len_next_nodes>1
            next_node_indx = randi(len_next_nodes,1);
            from = next_nodes(next_node_indx);
        else
            from = next_nodes;
        end
        k=k+1;
    end
    
    node_temp(k)=from;
    len_path = sum(~isnan(node_temp));
    node = node_temp(1:len_path);
    
    switch len_path
        case 1
            
            gibbs_block1_sub
        case 2
            
            gibbs_block2_sub
        case 3
            
            gibbs_block3_sub
        case 4
            
            gibbs_block4_sub
    end
    
    
    
    %Save dag and LL
    if toc>next_sample_time
        LL_retain(sampind)=LL;
        sampled_graphs{sampind} = sparse(logical(dag));
        
        sampind=sampind+1;
        next_sample_time = sampind*sample_duration;
        if next_sample_time >runt
            return;
        end
    end
    
    %         temp(t)=LL-score_dag(dag,LS,pa_limit,dummy,padim);
    %         if abs(temp(t))>0.5
    %             disp(temp(t))
    %             keyboard
    %         end
    
end


%%%%%%%%%


function [pa_loc,postPend1,n1_ls]=drawpa3(pb_H0,n1_Spa_kj1_kj2,postPend1,u)
% search move destination specified by random number u

n1_ls=[];
pa_loc=[];
postP=postPend1+cumsum(pb_H0);
if isempty(postP)
    return;
end
if  postP(end)<u
    postPend1=postP(end);
else
    if  postP(1)>u
        pa_loc = 1;
    else
        pa_loc = find(postP <u,1,'last')+1;
        postPend1=postP(pa_loc-1);
    end
    
    n1_ls = n1_Spa_kj1_kj2(pa_loc);
    
end



function pa_loc=drawpa3_1(pb_H0,postPend1,u)
% search move destination specified by random number u

postP=postPend1+cumsum(pb_H0);

if  postP(1)>u
    pa_loc = 1;
else
    pa_loc = find(postP <u,1,'last')+1;
    postPend1=postP(pa_loc-1);
end

