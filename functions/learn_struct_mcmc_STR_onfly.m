function [sampled_graphs, LL_retain,accept_ratio] = learn_struct_mcmc_STR_onfly(varargin)
%learn_struct_mcmc_STR learning network structure by mcmc
% INPUTS:
%	nsamples:	the number of samples collected from the chain
%   init_dag: initail network structure
%   LS: local score object
%   hprior: matrix of acceptable edges: 1-allowed, 0-banned
%   pa_limit: the number of maximum parents
%   runtime: running time
%   sample_duration: time intervel collecting samples
% OUTPUTS:
%	sampled_graphs:	sampled networks
%	LL_retain:	log likelihood of sampled networks
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License



args = varargin;
nargs = length(args);

Hprior=[];
runt=inf;
for i=1:2:nargs
    switch args{i},
        case 'nsamples',   nsamples = args{i+1}; %nsamples=100
        case 'init_dag',   dag = args{i+1};
        case 'pa_limit',    pa_limit = args{i+1};
        case 'hprior', Hprior = args{i+1};
        case 'runtime', runt=args{i+1};
        case 'sample_duration', sample_duration=args{i+1};
        case 'expdata', expdata=args{i+1};
    end
end

sampled_graphs = cell(1, nsamples);
dag = double(dag);


% reachability graph
n=length(dag);
[ai, bi] = find(dag==1);
cmat = zeros(n);
cmat(diag(ones(n,1))==1)=1;
for i=1:length(ai)
    cmat=cmat+cmat(:,ai(i))*cmat(bi(i),:);
end

% hard prior
if isempty(Hprior)
    Hprior = ones(n);
    Hprior = logical(setdiag(Hprior,0));
end


accept_ratio = zeros(1, nsamples);

LL_retain = -inf(nsamples,1);

numParent = sum(dag);
numEdge = sum(numParent);

LL = score_dag_onfly(dag,expdata); %+sum(log(pmat(dag==1)))+sum(log(1-pmat(dag==0)));

%%
[ops, nodes] = mk_nbrs_of_digraph(dag,cmat,pa_limit,Hprior,numEdge,numParent);

sampind=1;
next_sample_time = sampind*sample_duration;
tic
while 1
    
    dagp = dag;
    cmatp = cmat;
    
    % pick next move
    gIndx = randi(length(ops),1);
    source=nodes(gIndx,1);
    dest=nodes(gIndx,2);
    
    % update count matrix
    switch ops(gIndx)
        case 0
            dagp(source,dest)=0;
            cmatp = cmatp - cmatp(:,source)*cmatp(dest,:);
            numParent(dest) = numParent(dest)-1;
            numEdge = numEdge-1;
            
        case 1
            dagp(source,dest)=0;
            dagp(dest,source)=1;
            cmatp = cmatp - cmatp(:,source)*cmatp(dest,:);
            cmatp = cmatp + cmatp(:,dest)*cmatp(source,:);
            
            numParent(dest) = numParent(dest)-1;
            numParent(source) = numParent(source)+1;
            
            
        case 2
            dagp(source,dest)=1;
            cmatp = cmatp + cmatp(:,source)*cmatp(dest,:);
            
            numParent(dest) = numParent(dest)+1;
            numEdge = numEdge+1;
            
    end
    
    % get neighbors of the next move
    [opsp, nodesp] = mk_nbrs_of_digraph(dagp,cmatp,pa_limit,Hprior,numEdge,numParent);
    
    % calc local score
    paset_trt=find(dagp(:, dest))';
    LLnew = CcalcLS_one(expdata,dest,paset_trt);
    
    paset_trt=find(dag(:, dest))';
    LLold = CcalcLS_one(expdata,dest,paset_trt);
    
    if ops(gIndx)==1  % must also multiply in the changes to i's family
        
        paset_trt=find(dagp(:, source))';
        LLnew = LLnew+CcalcLS_one(expdata,source,paset_trt);
        
        tpaset_trt=find(dag(:, source))';
        LLold = LLold+CcalcLS_one(expdata,source,tpaset_trt);
    end
    
    bfactor = exp(LLnew - LLold);
    
    
    % decision
    R = bfactor * (length(ops) / length(opsp));
    u = rand(1,1);
    if u > R % reject the move
        switch ops(gIndx)
            case 0
                numParent(dest) = numParent(dest)+1;
                numEdge = numEdge+1;
                
            case 1
                numParent(dest) = numParent(dest)+1;
                numParent(source) = numParent(source)-1;
                
                
            case 2
                numParent(dest) = numParent(dest)-1;
                numEdge = numEdge-1;
                
        end
    else
        
        dag = dagp;
        ops = opsp;
        nodes = nodesp;
        cmat = cmatp;
        LL=LL+LLnew - LLold;
        
    end
    
    % save dag & LL
    if toc>next_sample_time
        LL_retain(sampind)=LL;
        sampled_graphs{sampind} = sparse(logical(dag));
        
        sampind=sampind+1;
        next_sample_time = sampind*sample_duration;
        if next_sample_time >runt
            return;
        end
    end
    
    
    %     if debug>1
    %         temp(t)=LL-score_dag4(dag,LS,pa_limit);
    %         if abs(temp(t))>0.5
    %             disp(temp(t))
    %         end
    %     end
    
end

end

