function [JTK_GRP_SIZE,JTK_NUM_GRPS,JTK_NUM_VALS,JTK_MAX,JTK_DIMS,JTK_SDV,JTK_EXV,JTK_EXACT,JTK_CP]=jtkdist(timepoints,reps)
% 
% clc;clear; close all;
% timepoints=912;
% reps=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];

len_reps=length(reps);
if(len_reps==timepoints)
    tim=reps;
end

normal=false;
len=length(tim);
JTK_GRP_SIZE=tim;
JTK_NUM_GRPS=len;
JTK_NUM_VALS=sum(tim);
JTK_MAX=(sum(tim)^2-sum(tim.^2))/2;
first_dim=sum(tim)*(sum(tim)-1)/2;
JTK_DIMS=[first_dim,1];

nn=JTK_NUM_VALS;
ns=JTK_GRP_SIZE;
%first_part=lfactorial(nn);
w=zeros();
for i=1:len
    w(i)=lfactorial(ns(i));
end
maxnlp=lfactorial(nn)-sum(w);
limit= log(realmax);

normal=normal |(maxnlp>(limit-1));

if(normal)
    JTK_VAR=(nn^2*(2*nn+3)-sum(ns.^2.*(2.*ns+3)))/72;
    JTK_SDV=sqrt(JTK_VAR);
    JTK_EXV=JTK_MAX/2;
    JTK_EXACT=false;
    JTK_CP=0;  %%no use in the future if normal is true
     return ;
end

MM=floor(JTK_MAX/2);
cf=num2cell(ones(1,MM+1));

sz=JTK_GRP_SIZE;
sz=sort(sz);         %%%sz is the same with size in R script
k=JTK_NUM_GRPS;
% N=s(k);
N=zeros();
if (k>2)
    for i=k:-1:2
        N(i-1)=sum(sz(i:k));  %%%%N is the same with N in R script
        
    end
end
N=N';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:k-1
   m=sz(i); %%m=0 or 1 
   n=N(i);
   
   if(n<MM)
       p=min(m+n,MM);  %%p<=n+1
       for t=(n+1):-1:p
           for u=(1+MM):-1:t
               cf{u}=cf{u}-cf{u-t+1};  %%In R script, cf[[u-t]], but u-t can be 0, so it may be wrong this part
           end
       end
   end
    
   Q=min(m,MM); %%m=0 or 1, so Q=0 or 1, should it be "Q=min(n,MM)"?
   for s=1:-1:Q
       for u=(1+s):MM
           cf{u}=cf{u}+cf{u-s};
       end
   end

end

cf=cell2mat(cf);  %transform cell into matrix, cf<-sapply(cf,sum)
cf_new=zeros();
part1=zeros();
part2=zeros();
part3=zeros();
part4=zeros();
if(mod(JTK_MAX,2))
    part2=[cf(MM:-1:1),0];
    part1=2*cf(MM+1)-part2;
    cf_new=[cf,part1];
else
    part4=[cf((MM-1):-1:1),0];
    part3=cf(MM+1)+cf(MM)-part4;
    cf_new=[cf,part3];
    
end

jtkcf=fliplr(cf_new); %reverse cf,  like rev() in R
len_jtkcf=length(jtkcf);
ajtkcf=(jtkcf(1:(len_jtkcf-1))+jtkcf(2:len_jtkcf))/2;

id=1+(0:(2*JTK_MAX));
cf_new=id;

len_id=length(id);
token1=1;
token2=1;
for ix=1:len_id
    if (mod(ix,2)) %%% "!!id%%2" is true
        cf_new(ix)=jtkcf(token1);       %odd value is jtkcf
        token1=token1+1;                 %
    else
        cf_new(ix)=ajtkcf(token2);      %even value is ajtkcf
        token2=token2+1;
    end
    
end


cp=cf_new/jtkcf(1);

JTK_CP=cp;
JTK_EXACT=TRUE;

 end









