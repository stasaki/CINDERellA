function [sampled_graphs, LL_retain] = learn_struct_mcmc_gibbs_block2_con(varargin)
%learn_struct_mcmc_gibbs_block2_con learning network structure by mcmc
% INPUTS:
%	nsamples:	the number of samples collected from the chain
%   init_dag: initail network structure
%   LS: local score object
%   hprior: matrix of acceptable edges: 1-allowed, 0-banned
%   pa_limit: the number of maximum parents
%   rev: try reversal edge move
%   pr: probablity to perform 2 parents block move
%   runtime: running time
%   sample_duration: time intervel collecting samples
% OUTPUTS:
%	sampled_graphs:	sampled networks
%	LL_retain:	log likelihood of sampled networks
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License



% debug=0;

args = varargin;
nargs = length(args);

% set default parameters
Pr = 0.25; % Probablity to performe 2 parent set blocking move
Hprior=[];
rev=0;
runt=inf;

% set paramters
for i=1:2:nargs
    switch args{i},
        case 'nsamples',   nsamples = args{i+1}; %nsamples=100
        case 'init_dag',   dag = args{i+1};
        case 'LS',       LS = args{i+1};
        case 'hprior', Hprior = args{i+1};
        case 'pa_limit',    pa_limit = args{i+1};
        case 'rev', rev=args{i+1};
        case 'pr', Pr = args{i+1};
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

% enumerate all node combination
nodecomb2=nchoosek(1:n,2);
nodecomb2=nodecomb2(randperm(size(nodecomb2,1)),:);
combl2=length(nodecomb2);
nperm=randperm(n);
bc2=0;
bc1=0;


temp2=repmat(dummy,1,pa_limit);

sampind=1;
next_sample_time = sampind*sample_duration;
tic
while 1
    
    %disp(t)
    
    
    
    if Pr > rand% modify parents for 2 node at the same time
        % select 2 nodes for modification
        if rev==1 % modify conneted nodes
            % store edge info
            [trt,src] = find(dag'&Hprior);
            nedge = length(src);
            if nedge==0 % if there is no edges in the network
                gibbs_block1_sub2 % perform one node modification
            else
                % pick one edge
                eIndx = randi(length(src),1); % if current graph is empty the code stops here
                node=[src(eIndx),trt(eIndx)];
                gibbs_block2_sub
            end
        else
            node=nodecomb2(mod(bc2,combl2)+1,:); %pick 2 nodes
            bc2=bc2+1;
            gibbs_block2_sub
        end
        
    else % modify parents for 1 node
        gibbs_block1_sub2
    end
    
    %Save dag and LL
    if toc>next_sample_time
        LL_retain(sampind)=LL;
        sampled_graphs{sampind} = sparse(logical(dag));
        
        sampind=sampind+1;
        next_sample_time = sampind*sample_duration;
        if toc >runt
            return;
        end
    end
    
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




