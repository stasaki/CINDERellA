%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License

%%


if flag_REV ==1 % if REV was performed previous step, need to update ops nodes.
    numParent = sum(dag);
    [ops, nodes] = mk_nbrs_of_digraph_del_add(dag,cmat,pa_limit,Hprior,numEdge,numParent);
    flag_REV = 0;
end




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
[opsp, nodesp] = mk_nbrs_of_digraph_del_add(dagp,cmatp,pa_limit,Hprior,numEdge,numParent);


% calc local score
patemp=LS{dest}.indx;
sctemp=LS{dest}.score;

paset_trt=find(dagp(:, dest))';
temp2(:)=dummy;
temp2(1:length(paset_trt))=paset_trt;
paset_trt=mysub2ind(padim,temp2);
LLnew = sctemp(ismembc2(paset_trt, patemp));



paset_trt=find(dag(:, dest))';
temp2(:)=dummy;
temp2(1:length(paset_trt))=paset_trt;
paset_trt=mysub2ind(padim,temp2);
LLold = sctemp(ismembc2(paset_trt, patemp));

if ops(gIndx)==1 

    patemp=LS{source}.indx;
    sctemp=LS{source}.score;
    
    paset_trt=find(dagp(:, source))';
    temp2(:)=dummy;
    temp2(1:length(paset_trt))=paset_trt;
    paset_trt=mysub2ind(padim,temp2);
    LLnew = LLnew+sctemp(ismembc2(paset_trt, patemp));
  
    tpaset_trt=find(dag(:, source))';
    temp2(:)=dummy;
    temp2(1:length(tpaset_trt))=tpaset_trt;
    paset_trt=mysub2ind(padim,temp2);
    LLold = LLold+sctemp(ismembc2(paset_trt, patemp));
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

