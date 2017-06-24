function array_one_column=get_low_tri(array)

sz=size(array);
N=sz(1);
array2=tril(array);

total_len=N*(N-1)/2;
array_one_column=zeros(total_len,1);
ix=0;

for i=1:(N-1)
%     num_eachcol=N-i;
    for j=(i+1):N
        ix=ix+1;
        array_one_column(ix,1)=array2(j,i);
    end
    
end
% ix


% for i=1:N
%     for j=1:N
%         if(array2(j,i)~=0)
%             ix=ix+1;
%             array_one_column(ix)=array2(j,i);
%             
%         else
%             continue;
%             
%         end
%     end
% end



end