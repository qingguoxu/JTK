function med=wilcoxon_median(x,mu)

x=x+mu;
len=length(x)*(1+length(x))/2;
y=zeros(len,1);
c=1;
for i=1:length(x)
    for j=i:length(x)
        y(c)=x(i)+x(j);
        c=c+1;
    end
end
y=sort(y)/2;


med=median(y);
end