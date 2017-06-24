function [p_output,s_output]=get_p_s(foosv,cgoosv,JTK_MAX,JTK_EXACT,JTK_EXV,JTK_SDV,JTK_CP)

        S=sum(foosv.*cgoosv);
        if(~S)
            p_output=1;
            s_output=0;
            
        else
            M=JTK_MAX;
            jtk=(abs(S)+M)/2;
            if(JTK_EXACT)
                jtki=1+2*jtk;
                p=2*JTK_CP(jtki);
                
            else
                temp=normcdf(-(jtk-1/2),-JTK_EXV,JTK_SDV);  %%%pnorm()
                p=2*temp;
                
            end
            p_output=p;
            s_output=S;
            
        end

end