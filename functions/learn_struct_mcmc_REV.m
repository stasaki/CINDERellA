function [sampled_graphs, LL_retain,accept_ratio] = learn_struct_mcmc_REV(varargin)
%learn_struct_mcmc_REV learning network structure by mcmc
% INPUTS:
%	nsamples:	the number of samples collected from the chain
%   init_dag: initail network structure
%   LS: local score object
%   hprior: matrix of acceptable edges: 1-allowed, 0-banned
%   pa_limit: the number of maximum parents
%   rev: probablity to perform reversal edge move
%   runtime: running time
%   sample_duration: time intervel collecting samples
% OUTPUTS:
%	sampled_graphs:	sampled networks
%	LL_retain:	log likelihood of sampled networks
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License



% debug = 0;

args = varargin;
nargs = length(args);

Hprior=[];
Pr = 1/15;
runt=inf;
for i=1:2:nargs
    switch args{i},
        case 'nsamples',   nsamples = args{i+1}; %nsamples=100
        case 'init_dag',   dag = args{i+1};
        case 'LS',       LS = args{i+1};
        case 'pa_limit',    pa_limit = args{i+1};
        case 'hprior', Hprior = args{i+1};
        case 'rev', Pr=args{i+1};
        case 'runtime', runt=args{i+1};
        case 'sample_duration', sample_duration=args{i+1};
        case 'sample_interval', sample_interval=args{i+1};
    end
end



% set number of iteration
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
num_accepts = 1;


flag_REV=1;
LL_retain = -inf(nsamples,1);
dummy=n+1;
padim=repmat(dummy,1,pa_limit);
nopavec=uint16(dummy:dummy+pa_limit-1);
temppa=nopavec;
numParent = sum(dag);
numEdge = sum(numParent);

LL = score_dag(dag,LS,pa_limit,dummy,padim); %+sum(log(pmat(dag==1)))+sum(log(1-pmat(dag==0)));

temp2=repmat(dummy,1,pa_limit);
%%
sampind=1;

if (exist('sample_duration')) %#ok<EXIST>
    next_sample_time = sampind*sample_duration;
end

stopByStep = false;
if (exist('sample_interval')) %#ok<EXIST>
    c_step=0;
    n_step=1;
    stopByStep = true;
end

tic
while 1
    
    %      disp(t)
    if numEdge==0 || Pr < rand
        MH_struct_sub
    else
        
        %% store edge info
        [trt,src] = find(dag'&Hprior); % find reversible edges
        nedge = length(src); % get the number of reversible edges
        if nedge==0 % if there's no reversible edges
            MH_struct_sub % do STR perfoming edege deletion or addtion
        else
            
            nedge_new = nedge;
            eIndx = randi(length(src),1); % select a edge to be reversed
            source=src(eIndx);
            target=trt(eIndx);
            
            MH_REVmove_sub
        end
    end
    
    if stopByStep
        c_step=c_step+1;
        
        if n_step==c_step
            LL_retain(sampind)=LL;
            sampled_graphs{sampind} = sparse(logical(dag));
            
            sampind=sampind+1;
            n_step = n_step+sample_interval;
            if sampind > nsamples
                return;
            end
        end
        continue;
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
    
    
    %     if debug>1
    %         temp(t)=LL-score_dag4(dag,LS,pa_limit);
    %         if abs(temp(t))>0.5
    %             disp(temp(t))
    %         end
    %     end
    
end

end

