function leap_year=isleapyear(year)
% year=2004;
sign=0;
if (rem(year,4)==0)
    sign=sign+1;
end
if (rem(year,100)==0)
    sign=sign-1;
end
if (rem(year,400)==0)
    sign=sign+1;
end
if sign==1
    
    leap_year=1; %%year is a leap year

else
    
    leap_year=0; %%year is not a leap year

end


end