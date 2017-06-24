function  output=jtk_rank(array)
%clear;clc;
%array=[32 ,43, 4, 1, 2, 3, 4,4,4,34,5,23,656,78,85,3,5,7,9];
len=length(array);
array2=sort(array);
tempoutput=zeros(1,len);
output=zeros(1,len);
token=1;

while token<=len
    
    endpoint=token;
    
    for i=(token+1):len
        
        if(array2(i)==array2(token))
            
            continue;
            
        else
            endpoint=i-1;
            break;
        end
    end
    
    if(endpoint~=token)
        tempoutput(token:endpoint)=mean(token:endpoint);
        token=endpoint+1;
    else
        
        tempoutput(token)=token;
        token=endpoint+1;
    end
    
    
    
end

for j=1:len
   for k=1:len
    if(array2(k)==array(j))
        
        output(j)=tempoutput(k);
    end
   end
    
end


%output

end