function [JTK_PERIOD,JTK_LAG,JTK_AMP]=jtkx(z,conf,JTK_DIMS,JTK_PERIODS,JTK_CGOOSV,JTK_SIGNCOS,JTK_MAX,JTK_EXACT,JTK_EXV,JTK_SDV,JTK_INTERVAL,JTK_AMPFACTOR,JTK_CP)
  
 JTK_PERIOD=0;JTK_LAG=0;JTK_AMP=0;
%  JTK_AMP_CI=0;
%  JTK_AMP_PVAL=0;
%  JTK_TAU=0;
 
 JTK_CJTK=jtkstat(z,JTK_DIMS,JTK_PERIODS,JTK_CGOOSV,JTK_MAX,JTK_EXACT,JTK_EXV,JTK_SDV,JTK_CP);
 cjtk=JTK_CJTK{1};
 pvals=cjtk(1,:);
 len=length(pvals);
 
 %%%%%%%%%%%%%%padj <- p.adjust(unlist(pvals),"bonf")%%
 padj=len*pvals;
 for i=1:len
     if(padj(i)>=1)
         padj(i)=1.0;
     end
 end
 %%%%%%%%%
 
 JTK_ADJP=min(padj);
 minpadj=min(padj);
 
 len_JTKADJP=length(JTK_ADJP);
 peris=zeros();
 inx=0;
 for i=1:len_JTKADJP
     if(JTK_ADJP(i)==minpadj)
         inx=inx+1;
         peris(inx)=i;
     end
 end
%  peris=1;
 pers=JTK_PERIODS(peris);
 
 lagis=zeros();
 inxx=0;
 len_padj=length(padj);
 for j=1:len_padj
     if(JTK_ADJP==padj(j))
         inxx=inxx+1;
         lagis(inxx)=j;
     else
         continue;
     end
 end
 
 len_lagis=length(lagis);
 temp_length=zeros();
 inxxx=0;
 for m=1:len_lagis
     inxxx=inxxx+1;
     temp_length(inxxx)=length(lagis(m));
 end
 
 count=sum(temp_length);
 
 bestper=0;
 bestlag=0;
 maxamp=0;
 maxamp_ci=zeros(2,1);
 maxamp_pval=0;
 
 for i=1:length(pers)
% i=1;
     per=pers(i);
     peri=peris(i);
     cjtk_jtkx=JTK_CJTK{peri};
     sc=JTK_SIGNCOS{peri};
     [nrow_sc,ncol_sc]=size(sc);
     temp_w=z(1:nrow_sc,1);
     median_temp_w=hlm(temp_w);
     w=(temp_w-median_temp_w)*JTK_AMPFACTOR;
     
     lagi=lagis(i);
     S_upper=cjtk_jtkx(2,lagi);
     s_little=sign(S_upper);
     if(~s_little)
         s_little=1;
     end
     
     temp_lag=(per+(1-s_little)*per/4-(lagi-1)/2);
     lag=mod(temp_lag,per);
     signcos=sc(:,lagi);
     tmp=s_little*w.*signcos;%%%
     
     
%      amp=median(tmp);
     amp=wilcoxon_median(tmp,0);
     if(amp>maxamp)
         maxamp=amp;
         bestper=per;
         bestlag=lag;
         
     end
     
     JTK_PERIOD=JTK_INTERVAL*bestper;
     JTK_LAG=JTK_INTERVAL*bestlag;
     JTK_AMP=max(0,maxamp);
     JTK_TAU=abs(S_upper)/JTK_MAX;
     
%      [pjkt hjtk]=ranksum(tmp,'alpha',0.1,'method','approximate');



 end
 
 
end