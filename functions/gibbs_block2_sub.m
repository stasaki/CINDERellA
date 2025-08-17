%% sub code for learn_struct_mcmc_gibbs_block12_v4
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License


pa_node1 = find(dag(:,node(1)))'; %get current parents for node1
pa_node2 = find(dag(:,node(2)))'; %get current parents for node2

%orphaning node1 and node2 and updating path count matrix(cmat)
for i=1:length(pa_node1)
    dag(pa_node1(i),node(1))=0;
    cmat = cmat - cmat(:,pa_node1(i))*cmat(node(1),:);
end
for i=1:length(pa_node2)
    dag(pa_node2(i),node(2))=0;
    cmat = cmat - cmat(:,pa_node2(i))*cmat(node(2),:);
end


n1_lookup=LS{node(1)}.lookup;
n2_lookup=LS{node(2)}.lookup;

%retrive parent candidates
pa_kj1c = find(cmat(node(1),:)~=0); %descendant for node1
pa_kj2c = find(cmat(node(2),:)~=0); %descendant for node2


%common descendant for node1 and node2
pa_kj1c_kj2c=unionKPM(pa_kj1c,pa_kj2c);

%compute parent sets
%find possible parent candidates TYPE0 for node1 and node2
n1_indx1=find(~any(n1_lookup(:,pa_kj1c_kj2c),2));
n2_indx1=find(~any(n2_lookup(:,pa_kj1c_kj2c),2));

%find possible parent candidates TYPE1 for node1
n1_indx2=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj2c),2));

%find possible parent candidates TYPE1 for node2 which
n2_indx2=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj1c),2));


%get local score for possible child and parents combination for node1
sctemp=LS{node(1)}.score;
n1_LS_Spa_kj1_kj2=sctemp(n1_indx1);
n1_LS_Spa_kj1_kj2c=sctemp(n1_indx2);


temp2(:)=dummy;
temp2(1:length(pa_node1))=pa_node1;
paset_trt=mysub2ind(padim,temp2);
n1_ls_old = sctemp(ismembc2(paset_trt, LS{node(1)}.indx));


%get local score for possible child and parents combination for node2
sctemp=LS{node(2)}.score;
n2_LS_Spa_kj1_kj2=sctemp(n2_indx1);
n2_LS_Spa_kj2_kj1c=sctemp(n2_indx2);

temp2(:)=dummy;
temp2(1:length(pa_node2))=pa_node2;
paset_trt=mysub2ind(padim,temp2);
n2_ls_old = sctemp(ismembc2(paset_trt, LS{node(2)}.indx));


logsumtemp=logsumexp(n2_LS_Spa_kj1_kj2,1);
pb_H0=n1_LS_Spa_kj1_kj2+logsumtemp;

if isempty(n1_LS_Spa_kj1_kj2c)
    pb_H1=-Inf;
else
    pb_H1=n1_LS_Spa_kj1_kj2c+logsumtemp;
end

if isempty(n2_LS_Spa_kj2_kj1c)
    pb_H2=-Inf+n1_LS_Spa_kj1_kj2;
else
    pb_H2=n1_LS_Spa_kj1_kj2+logsumexp(n2_LS_Spa_kj2_kj1c,1);
end
normconst=logsumexp([pb_H0;pb_H1;pb_H2],1);


u=rand;
postPend=0;% store n1_ls and n2_ls for debug use
[n1_loc,postPend,n1_ls]=drawpa3(exp(pb_H0-normconst),n1_LS_Spa_kj1_kj2,postPend,u);
if isempty(n1_loc)
    [n1_loc,postPend,n1_ls]=drawpa3(exp(pb_H1-normconst),n1_LS_Spa_kj1_kj2c,postPend,u);
    if isempty(n1_loc)
        [n1_loc,postPend,n1_ls]=drawpa3(exp(pb_H2-normconst),n1_LS_Spa_kj1_kj2,postPend,u);
        n1_pa = LS{node(1)}.parent(n1_indx1(n1_loc),:);
        
        n2_loc=drawpa3_1(exp(n1_ls+n2_LS_Spa_kj2_kj1c-normconst),postPend,u);
        n2_pa = LS{node(2)}.parent(n2_indx2(n2_loc),:);
        n2_ls = n2_LS_Spa_kj2_kj1c(n2_loc);
    else
        n1_pa = LS{node(1)}.parent(n1_indx2(n1_loc),:);
        
        n2_loc=drawpa3_1(exp(n1_ls+n2_LS_Spa_kj1_kj2-normconst),postPend,u);
        n2_pa = LS{node(2)}.parent(n2_indx1(n2_loc),:);
        n2_ls = n2_LS_Spa_kj1_kj2(n2_loc);
    end
else
    n1_pa = LS{node(1)}.parent(n1_indx1(n1_loc),:);
    n2_loc=drawpa3_1(exp(n1_ls+n2_LS_Spa_kj1_kj2-normconst),postPend,u);
    n2_pa = LS{node(2)}.parent(n2_indx1(n2_loc),:);
    n2_ls = n2_LS_Spa_kj1_kj2(n2_loc);
end



%Update dag and cmat
n1_pa=rmzero(n1_pa,dummy);
n2_pa=rmzero(n2_pa,dummy);
dag(n1_pa,node(1))=1;
dag(n2_pa,node(2))=1;


for j=1:length(n1_pa)
    cmat = cmat + cmat(:,n1_pa(j))*cmat(node(1),:);
end
for j=1:length(n2_pa)
    cmat = cmat + cmat(:,n2_pa(j))*cmat(node(2),:);
end

%Update current LL
LL = LL + n1_ls + n2_ls - n1_ls_old - n2_ls_old;

