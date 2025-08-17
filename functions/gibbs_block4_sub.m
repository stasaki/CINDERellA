%% sub code for learn_struct_mcmc_gibbs_block12_v2
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License



pa_node1 = find(dag(:,node(1)))'; %get current parents for node1
pa_node2 = find(dag(:,node(2)))'; %get current parents for node2
pa_node3 = find(dag(:,node(3)))'; %get current parents for node3
pa_node4 = find(dag(:,node(4)))'; %get current parents for node4

%%orphaning node1 and node2 and updating path count matrix(cmat)
for i=1:length(pa_node1)
    dag(pa_node1(i),node(1))=0;
    cmat = cmat - cmat(:,pa_node1(i))*cmat(node(1),:);
end
for i=1:length(pa_node2)
    dag(pa_node2(i),node(2))=0;
    cmat = cmat - cmat(:,pa_node2(i))*cmat(node(2),:);
end
for i=1:length(pa_node3)
    dag(pa_node3(i),node(3))=0;
    cmat = cmat - cmat(:,pa_node3(i))*cmat(node(3),:);
end
for i=1:length(pa_node4)
    dag(pa_node4(i),node(4))=0;
    cmat = cmat - cmat(:,pa_node4(i))*cmat(node(4),:);
end
%%
n1_lookup=LS{node(1)}.lookup;
n2_lookup=LS{node(2)}.lookup;
n3_lookup=LS{node(3)}.lookup;
n4_lookup=LS{node(4)}.lookup;

%retrive parent candidates
pa_kj1c = find(cmat(node(1),:)~=0); %descendant for node1
pa_kj2c = find(cmat(node(2),:)~=0); %descendant for node2
pa_kj3c = find(cmat(node(3),:)~=0); %descendant for node3
pa_kj4c = find(cmat(node(4),:)~=0); %descendant for node4

%union descendant for node1 and node2 and node3
pa_kj1c_kj2c_kj3c=unionKPM(pa_kj4c,unionKPM(pa_kj3c,unionKPM(pa_kj1c,pa_kj2c)));%back here later

%%compute parent sets
%% find possible parent candidates TYPE0 for node1, node2 and node3
n1_4indx{1}=find(~any(n1_lookup(:,pa_kj1c_kj2c_kj3c),2));
n2_4indx{1}=find(~any(n2_lookup(:,pa_kj1c_kj2c_kj3c),2));
n3_4indx{1}=find(~any(n3_lookup(:,pa_kj1c_kj2c_kj3c),2));
n4_4indx{1}=find(~any(n4_lookup(:,pa_kj1c_kj2c_kj3c),2));

%% find possible parent candidates TYPE1 for node1
%node1<-node2
n1_4indx{2}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj2c),2)&(~any(n1_lookup(:,pa_kj3c),2))&(~any(n1_lookup(:,pa_kj4c),2)));
%node1<-node3
n1_4indx{3}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj3c),2)&(~any(n1_lookup(:,pa_kj2c),2))&(~any(n1_lookup(:,pa_kj4c),2)));
%node1<-node4
n1_4indx{4}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj4c),2)&(~any(n1_lookup(:,pa_kj2c),2))&(~any(n1_lookup(:,pa_kj3c),2)));

%find possible parent candidates TYPE1 for node2
%node2<-node1
n2_4indx{2}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj1c),2)&(~any(n2_lookup(:,pa_kj3c),2))&(~any(n2_lookup(:,pa_kj4c),2)));
%node2<-node3
n2_4indx{3}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj3c),2)&(~any(n2_lookup(:,pa_kj1c),2))&(~any(n2_lookup(:,pa_kj4c),2)));
%node2<-node4
n2_4indx{4}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj4c),2)&(~any(n2_lookup(:,pa_kj1c),2))&(~any(n2_lookup(:,pa_kj3c),2)));


%find possible parent candidates TYPE1 for node3
%node3<-node1
n3_4indx{2}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj1c),2)&(~any(n3_lookup(:,pa_kj2c),2))&(~any(n3_lookup(:,pa_kj4c),2)));
%node3<-node2
n3_4indx{3}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj2c),2)&(~any(n3_lookup(:,pa_kj1c),2))&(~any(n3_lookup(:,pa_kj4c),2)));
%node3<-node4
n3_4indx{4}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj4c),2)&(~any(n3_lookup(:,pa_kj1c),2))&(~any(n3_lookup(:,pa_kj2c),2)));

%find possible parent candidates TYPE1 for node4
%node4<-node1
n4_4indx{2}=find((~any(n4_lookup(:,pa_kj4c),2))&any(n4_lookup(:,pa_kj1c),2)&(~any(n4_lookup(:,pa_kj2c),2))&(~any(n4_lookup(:,pa_kj3c),2)));
%node4<-node2
n4_4indx{3}=find((~any(n4_lookup(:,pa_kj4c),2))&any(n4_lookup(:,pa_kj2c),2)&(~any(n4_lookup(:,pa_kj1c),2))&(~any(n4_lookup(:,pa_kj3c),2)));
%node4<-node3
n4_4indx{4}=find((~any(n4_lookup(:,pa_kj4c),2))&any(n4_lookup(:,pa_kj3c),2)&(~any(n4_lookup(:,pa_kj1c),2))&(~any(n4_lookup(:,pa_kj2c),2)));

%% find possible parent candidates TYPE2 for node1
%node1<-(node2,node3)
n1_4indx{5}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj2c),2)&(any(n1_lookup(:,pa_kj3c),2))&(~any(n1_lookup(:,pa_kj4c),2)));
%node1<-(node2,node4)
n1_4indx{6}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj2c),2)&(any(n1_lookup(:,pa_kj4c),2))&(~any(n1_lookup(:,pa_kj3c),2)));
%node1<-(node3,node4)
n1_4indx{7}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj3c),2)&(any(n1_lookup(:,pa_kj4c),2))&(~any(n1_lookup(:,pa_kj2c),2)));


%find possible parent candidates TYPE2 for node2
%node2<-(node1,node3)
n2_4indx{5}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj1c),2)&(any(n2_lookup(:,pa_kj3c),2))&(~any(n2_lookup(:,pa_kj4c),2)));
%node2<-(node1,node4)
n2_4indx{6}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj1c),2)&(any(n2_lookup(:,pa_kj4c),2))&(~any(n2_lookup(:,pa_kj3c),2)));
%node2<-(node3,node4)
n2_4indx{7}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj3c),2)&(any(n2_lookup(:,pa_kj4c),2))&(~any(n2_lookup(:,pa_kj1c),2)));


%find possible parent candidates TYPE2 for node3
%node3<-(node1,node2)
n3_4indx{5}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj1c),2)&(any(n3_lookup(:,pa_kj2c),2))&(~any(n3_lookup(:,pa_kj4c),2)));
%node3<-(node1,node4)
n3_4indx{6}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj1c),2)&(any(n3_lookup(:,pa_kj4c),2))&(~any(n3_lookup(:,pa_kj2c),2)));
%node3<-(node2,node4)
n3_4indx{7}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj4c),2)&(any(n3_lookup(:,pa_kj4c),2))&(~any(n3_lookup(:,pa_kj1c),2)));


%find possible parent candidates TYPE2 for node4
%node4<-(node1,node2)
n4_4indx{5}=find((~any(n4_lookup(:,pa_kj4c),2))&any(n4_lookup(:,pa_kj1c),2)&(any(n4_lookup(:,pa_kj2c),2))&(~any(n4_lookup(:,pa_kj3c),2)));
%node4<-(node1,node3)
n4_4indx{6}=find((~any(n4_lookup(:,pa_kj4c),2))&any(n4_lookup(:,pa_kj1c),2)&(any(n4_lookup(:,pa_kj3c),2))&(~any(n4_lookup(:,pa_kj2c),2)));
%node4<-(node2,node3)
n4_4indx{7}=find((~any(n4_lookup(:,pa_kj4c),2))&any(n4_lookup(:,pa_kj2c),2)&(any(n4_lookup(:,pa_kj3c),2))&(~any(n4_lookup(:,pa_kj1c),2)));

%% find possible parent candidates TYPE3 for node1
%node1<-(node2,node3,node4)
n1_4indx{8}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj2c),2)&(any(n1_lookup(:,pa_kj3c),2))&(any(n1_lookup(:,pa_kj4c),2)));

%find possible parent candidates TYPE2 for node2
%node2<-(node1,node3,node4)
n2_4indx{8}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj1c),2)&(any(n2_lookup(:,pa_kj3c),2))&(any(n2_lookup(:,pa_kj4c),2)));

%find possible parent candidates TYPE2 for node3
%node3<-(node1,node2,node4)
n3_4indx{8}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj1c),2)&(any(n3_lookup(:,pa_kj2c),2))&(any(n3_lookup(:,pa_kj4c),2)));

%find possible parent candidates TYPE2 for node4
%node4<-(node1,node2,node3)
n4_4indx{8}=find((~any(n4_lookup(:,pa_kj4c),2))&any(n4_lookup(:,pa_kj1c),2)&(any(n4_lookup(:,pa_kj2c),2))&(any(n4_lookup(:,pa_kj3c),2)));

%% get local score for possible child and parents combination for node1

sctemp=LS{node(1)}.score;
for i=1:8
    n1_4LS_indx{i}=sctemp(n1_4indx{i});
end

temp2(:)=dummy;
temp2(1:length(pa_node1))=pa_node1;
paset_trt=mysub2ind(padim,temp2);
n1_ls_old = sctemp(ismembc2(paset_trt, LS{node(1)}.indx));


%% get local score for possible child and parents combination for node2

sctemp=LS{node(2)}.score;
for i=1:8
    n2_4LS_indx{i}=sctemp(n2_4indx{i});
end

temp2(:)=dummy;
temp2(1:length(pa_node2))=pa_node2;
paset_trt=mysub2ind(padim,temp2);
n2_ls_old = sctemp(ismembc2(paset_trt, LS{node(2)}.indx));


%% get local score for possible child and parents combination for node3

sctemp=LS{node(3)}.score;
for i=1:8
    n3_4LS_indx{i}=sctemp(n3_4indx{i});
end

temp2(:)=dummy;
temp2(1:length(pa_node3))=pa_node3;
paset_trt=mysub2ind(padim,temp2);
n3_ls_old = sctemp(ismembc2(paset_trt, LS{node(3)}.indx));


%% get local score for possible child and parents combination for node4

sctemp=LS{node(4)}.score;
for i=1:8
    n4_4LS_indx{i}=sctemp(n4_4indx{i});
end


temp2(:)=dummy;
temp2(1:length(pa_node4))=pa_node4;
paset_trt=mysub2ind(padim,temp2);
n4_ls_old = sctemp(ismembc2(paset_trt, LS{node(4)}.indx));


%% calc normalizing constant


for i=1:48
    pb_4n3{i}=n3_4LS_indx{mar4code3(i,1)}+logsumexp(n4_4LS_indx{mar4code3(i,2)},1);
end


for i=1:200
    pb_4n2{i}=n2_4LS_indx{mar4code2(i,1)}+logsumexp(pb_4n3{mar4code2(i,2)},1);
end


for i=1:543
    pb_4n1{i}=n1_4LS_indx{mar4code1(i,1)}+logsumexp(pb_4n2{mar4code1(i,2)},1);
end

normconst=logsumexp(cat(1,pb_4n1{:}),1);

n4_loc=[];
while isempty(n4_loc)
    %Draw new parent
    u=rand;
    postPend=0;% store n1_ls and n2_ls for debug use
    
    % draw parant for n1
    for i=1:543
        [n1_loc,postPend,n1_ls]=drawpa3(exp(pb_4n1{i}-normconst),n1_4LS_indx{mar4code1(i,1)},postPend,u);
        
        if ~isempty(n1_loc)
            
            sp_loc=i;
            n1_pa = LS{node(1)}.parent(n1_4indx{mar4code1(i,1)}(n1_loc),:);
            break;
        end
    end
    [n2_loc,postPend,n2_ls]=drawpa3(exp(n1_ls+pb_4n2{mar4code(sp_loc,2)}-normconst),n2_4LS_indx{net4code(sp_loc,2)},postPend,u);
    [n3_loc,postPend,n3_ls]=drawpa3(exp(n1_ls+n2_ls+pb_4n3{mar4code(sp_loc,3)}-normconst),n3_4LS_indx{net4code(sp_loc,3)},postPend,u);
    [n4_loc,~,n4_ls]=drawpa3(exp(n1_ls+n2_ls+n3_ls+n4_4LS_indx{mar4code(sp_loc,4)}-normconst),n4_4LS_indx{net4code(sp_loc,4)},postPend,u);
end

n2_pa = LS{node(2)}.parent(n2_4indx{net4code(sp_loc,2)}(n2_loc),:);
n3_pa = LS{node(3)}.parent(n3_4indx{net4code(sp_loc,3)}(n3_loc),:);
n4_pa = LS{node(4)}.parent(n4_4indx{net4code(sp_loc,4)}(n4_loc),:);

%Update dag and cmat
n1_pa=rmzero(n1_pa,dummy);
n2_pa=rmzero(n2_pa,dummy);
n3_pa=rmzero(n3_pa,dummy);
n4_pa=rmzero(n4_pa,dummy);
dag(n1_pa,node(1))=1;
dag(n2_pa,node(2))=1;
dag(n3_pa,node(3))=1;
dag(n4_pa,node(4))=1;


for j=1:length(n1_pa)
    cmat = cmat + cmat(:,n1_pa(j))*cmat(node(1),:);
end
for j=1:length(n2_pa)
    cmat = cmat + cmat(:,n2_pa(j))*cmat(node(2),:);
end
for j=1:length(n3_pa)
    cmat = cmat + cmat(:,n3_pa(j))*cmat(node(3),:);
end
for j=1:length(n4_pa)
    cmat = cmat + cmat(:,n4_pa(j))*cmat(node(4),:);
end


%Update current LL
LL = LL + n1_ls + n2_ls + n3_ls + n4_ls- n1_ls_old - n2_ls_old - n3_ls_old -n4_ls_old;

