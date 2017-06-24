function output_median=hlm(w)

zz=outer_add(w,w);
sz=size(zz);
N=sz(1);
zz2=tril(zz,0);

total_len=N*(N+1)/2;
array_one_column=zeros(total_len,1);
ix=0;

for i=1:N  %%column
%     num_eachcol=N-i+1;
    for j=i:N
        ix=ix+1;
        array_one_column(ix,1)=zz2(j,i);
    end
    
end

temp_median=median(array_one_column);

output_median=temp_median/2;

end