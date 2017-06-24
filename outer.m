function arrayAB=outer(arrayA,arrayB)

if(length(arrayA)~=length(arrayB))
    display('errors in outer');
    return;
    
else
    len=length(arrayA);
    arrayAB=zeros(len,len);
    
    for i=1:len
        for j=1:len
            arrayAB(i,j)=arrayA(i)-arrayB(j);
        end
    end
    
end

end