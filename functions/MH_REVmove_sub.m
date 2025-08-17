%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License


%% store parents of target
paset_trt = find(dag(:,target))';
paset_src = find(dag(:,source))';


src_lookup=LS{source}.lookup;
src_score4=LS{source}.score;
trt_lookup=LS{target}.lookup;
trt_score4=LS{target}.score;

%% remove parents from source and target
dagp = dag;


dagp(:,source)=0;
dagp(:,target)=0;

% update cmat
cmatp = cmat;

for i=1:length(paset_src)
    nedge_new=nedge_new-1;
    cmatp = cmatp - cmatp(:,paset_src(i))*cmatp(source,:);
end

% store descendants of node source in dagpp
ds_src_dagpp = (cmatp(source,:)>0);

for i=1:length(paset_trt)
    nedge_new=nedge_new-1;
    cmatp = cmatp - cmatp(:,paset_trt(i))*cmatp(target,:);
end



%% store all descendants of source and target

ds_src = (cmatp(source,:)>0);%Hprior really needed?
ds_trt = (cmatp(target,:)>0);

%% pick and set new parent set for source

%         determin candidates of parent sets


pcIndx=find((~any(src_lookup(:,ds_src),2))&any(src_lookup(:,target),2));
postLL_pasets_src = src_score4(pcIndx);



%Draw new parent
u=rand;
partf1 = postLL_pasets_src;
expLL=exp(postLL_pasets_src-max(postLL_pasets_src));

postP=cumsum(expLL./sum(expLL));
if  postP(1)>u
    pa_loc = 1;
else
    pa_loc = find(postP <u,1,'last')+1;
end

LLnew = postLL_pasets_src(pa_loc);% for efficiency check, dont have to calc
% set new parent for source
new_pa = rmzero(LS{source}.parent(pcIndx(pa_loc),:),dummy);
dagp(new_pa,source)=1;

% update cmat
for i=1:length(new_pa)
    nedge_new=nedge_new+1;
    cmatp = cmatp + cmatp(:,new_pa(i))*cmatp(source,:);
end

%% store all descendatns of target
ds_trt_new = (cmatp(target,:)>0);

%% pick and set new parent set for target


%determin candidates of parent sets
pcIndx=find(~any(trt_lookup(:,ds_trt_new),2));
postLL_pasets_src = trt_score4(pcIndx);


%Draw new parent
u=rand;
partf2 = postLL_pasets_src;
expLL = exp(postLL_pasets_src-max(postLL_pasets_src));

postP=cumsum(expLL./sum(expLL));
if  postP(1)>u
    pa_loc = 1;
else
    pa_loc = find(postP <u,1,'last')+1;
end

LLnew = LLnew+postLL_pasets_src(pa_loc);% for efficiency check, dont have to calc
% set new parent for source
new_pa = rmzero(LS{target}.parent(pcIndx(pa_loc),:),dummy);
dagp(new_pa,target)=1;

% update cmat
for i=1:length(new_pa)
    nedge_new=nedge_new+1;
    cmatp = cmatp + cmatp(:,new_pa(i))*cmatp(target,:);
end


%% compute acceptance probability


pcIndx=(~any(trt_lookup(:,ds_trt),2))&any(trt_lookup(:,source),2);
partf3 = trt_score4(pcIndx);


pcIndx=find(~any(src_lookup(:,ds_src_dagpp),2));
partf4 = src_score4(pcIndx);


max_f1_f4=max([partf1;partf4]);
max_f2_f3=max([partf2;partf3]);
f1 = sum(exp(partf1-max_f1_f4));
f4 = sum(exp(partf4-max_f1_f4));
f2 = sum(exp(partf2-max_f2_f3));
f3 = sum(exp(partf3-max_f2_f3));


R=nedge_new*f1*f2/nedge/f3/f4;

%% decision

u = rand(1,1);
if u > R % reject the move

else
    
    dag = dagp;
    cmat = cmatp;
    
    numEdge = length(find(dag));

    
    temp2(:)=dummy;
    temp2(1:length(paset_trt))=paset_trt;
    paset_trt=mysub2ind(padim,temp2);
    LLold = trt_score4(ismembc2(paset_trt, LS{target}.indx));

    temp2(:)=dummy;
    temp2(1:length(paset_src))=paset_src;
    paset_trt=mysub2ind(padim,temp2);
    LLold = LLold+src_score4(ismembc2(paset_trt, LS{source}.indx));
    
    
    LL=LL+LLnew - LLold;

    num_accepts = num_accepts + 1;
    flag_REV = 1; % flag to indicate REV was performed
end