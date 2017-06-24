function [JTK_CGOOSV,JTK_SIGNCOS,JTK_PERIODS,JTK_INTERVAL]=jtk_init(periods,interval,JTK_GRP_SIZE,JTK_NUM_GRPS,JTK_NUM_VALS,JTK_DIMS,JTK_PIHAT)

%%%%%comment part
% clc;clear;close all;
% reps=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% tim=reps;
% periods=144:144;
% interval=0.16667;
% JTK_GRP_SIZE=reps;
% JTK_NUM_GRPS=length(reps);
% JTK_NUM_VALS=sum(reps);
% first_dim=sum(tim)*(sum(tim)-1)/2;
% JTK_DIMS=[first_dim,1];
% JTK_PIHAT= 3.1416;
%%%%%%

JTK_INTERVAL=interval;
JTK_PERIODS=periods;
JTK_PERFACTOR=ones(1,periods);

tim=JTK_GRP_SIZE;
timepoints=JTK_NUM_GRPS;
timerange=0:(timepoints-1);
len=length(tim);

JTK_CGOOSV=cell(1);
JTK_SIGNCOS=cell(1);
 
for i=1:length(periods)
   
    period=periods(i);
    ncol=period;    %%144
    nrow=JTK_DIMS(1); %%391170
    tmp_jtk_cgoosv=zeros(nrow,ncol);
    
    
    time2angle=2*JTK_PIHAT/period;
    theta=timerange*time2angle;
    cos_v=cos(theta);
    cos_r_temp=jtk_rank(cos_v);
    
%%%%%%%%this is for 'cos.r <- rep(cos.r,ti=tim)'  

    cos_r=zeros(JTK_NUM_VALS,1);
    ix=0;
    for j=1:len
        ti=tim(j);
        if(ti~=0)
            ix=ix+1;
            cos_r(ix)=cos_r_temp(j);
            
        else
            continue;
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    output_outer=outer(cos_r,cos_r);
    cgoos_first=sign(output_outer);
    cgoos=get_low_tri(cgoos_first);
    cgoosv=reshape(cgoos,JTK_DIMS(1),JTK_DIMS(2));
    
   
  %  JTK_CGOOSV{i}(:,1)=zeros(nrow,ncol);
  %  JTK_CGOOSV{i}=cgoosv;
  
    tmp_jtk_cgoosv(:,1)=cgoosv; 
    
    
    cycles=floor(timepoints/period);
    range=1:(cycles*period);
    cos_s_temp=sign(cos_v(range));
    %%%%%%%%this is for 'cos.s <- rep(cos.s,ti=tim[range])'  
    
    tim2=tim(range);
    JTK_NUM_VALS_2=sum(tim2);
    len2=length(tim2);
    cos_s=zeros(JTK_NUM_VALS_2,1);
    ixx=0;
    for j=1:len2
        ti2=tim2(j);
        if(ti2~=0)
            ixx=ixx+1;
            cos_s(ixx)=cos_s_temp(j);
            
        else
            continue;
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    nrow2=length(cos_s);
    ncol2=period;
    tmp_jtk_signcos=zeros(nrow2,ncol2);
    tmp_jtk_signcos(:,1)=cos_s; 
    
    for m=2:period
        delta.theta=(m-1)*time2angle/2;
        cos_v2=cos(theta+delta.theta);
        cos_r_temp2=jtk_rank(cos_v2);
        %%%%%%%%this is for 'cos.r <- rep(cos.r,ti=tim)'  

        cos_r2=zeros(JTK_NUM_VALS,1);
        inx=0;
        for n=1:len
            ti3=tim(n);
            if(ti3~=0)
                inx=inx+1;
                cos_r2(inx)=cos_r_temp2(n);
            
            else
                continue;
            end
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        output_outer2=outer(cos_r2,cos_r2);
        cgoos_first2=sign(output_outer2);
        cgoos2=get_low_tri(cgoos_first2);
        cgoosv2=reshape(cgoos2,JTK_DIMS(1),JTK_DIMS(2));
        tmp_jtk_cgoosv(:,m)=cgoosv2;
        
        cos_s_temp2=sign(cos_v2(range));
        %%%%%%%%this is for 'cos.s <- rep(cos.s,ti=tim[range])'  
        len2=length(tim2);
        cos_s2=zeros(JTK_NUM_VALS_2,1);
        inxx=0;
        for k=1:len2
            ti=tim2(k);
            if(ti~=0)
                inxx=inxx+1;
                cos_s2(inxx)=cos_s_temp2(k);
            
            else
                continue;
            end
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        tmp_jtk_signcos(:,m)=cos_s2;
    end
    JTK_CGOOSV{1}(:,:)= tmp_jtk_cgoosv(:,:);
    JTK_SIGNCOS{1}(:,:)= tmp_jtk_signcos(:,:);
end



 end

