function JTK_CJTK=jtkstat(z,JTK_DIMS,JTK_PERIODS,JTK_CGOOSV,JTK_MAX,JTK_EXACT,JTK_EXV,JTK_SDV,JTK_CP)

temp_output=outer(z,z);
foosv_temp1=sign(temp_output);
foosv_temp2=get_low_tri(foosv_temp1);
foosv=reshape(foosv_temp2,JTK_DIMS);   %%dim(foosv) <- JTK.DIMS

JTK_CJTK=cell(1);
temp_ps_array=zeros(2,144);
len=length(JTK_PERIODS);
for i=1:len
    jtkcgoosv=JTK_CGOOSV{i};
    [nrow,ncol]=size(jtkcgoosv);
    
    for j=1:ncol
        cg=jtkcgoosv(:,j);
        [p,S]=get_p_s(foosv,cg,JTK_MAX,JTK_EXACT,JTK_EXV,JTK_SDV,JTK_CP);
        temp_ps_array(1,j)=p;
        temp_ps_array(2,j)=S;  %%capital S
    end
    JTK_CJTK{i}=temp_ps_array;
    
    
end


end %%end for jtkstat function
