function S2=s2hss(S2t,IA2,dis,wr,n,x11,x00,y11,y00)
[mm,nn]=size(dis);
S2=[];
 s1=zeros(nn,1);
for r=1:nn
    for xi=1:n
        for yi=1:n
            if IA2(xi,yi)==1
                minwq=min(abs(xi-x00),(n-abs(xi-x00)));
                minhq=min(abs(yi-y00),(n-abs(yi-y00)));  
                rrq=sqrt(minwq^2+minhq^2);
                minwp=min(abs(xi-x11),(n-abs(xi-x11)));
                minhp=min(abs(yi-y11),(n-abs(yi-y11)));                        
                rrp=sqrt(minwp^2+minhp^2);
                if (rrq==dis(r))&(rrp~=dis(r))
                    s1(r)=s1(r)+1;
                elseif (rrq~=dis(r))&(rrp==dis(r))
                    s1(r)=s1(r)-1;
                else
                    s1(r)=s1(r)+0;
                end
            end
        end
    end
    if dis(r)==0
        s1(r)=s1(r)+1;
    end     
      minwpq=min(abs(x00-x11),(n-abs(x00-x11)));
      minhpq=min(abs(y00-y11),(n-abs(y00-y11)));                        
      rrpq=sqrt(minwpq^2+minhpq^2);
     if rrpq==dis(r)
         s1(r)=s1(r)-1;
     end                        
    S2(r)=S2t(r)+2*s1(r)/(n*n*wr(r));%得到相关函数。
end  
end
