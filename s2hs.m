function S2=s2hs(dis,IA1,wr,n)
%dis为距离r矩阵
%IA1为要求相关函数的结构的表示矩阵
%wr表示距离一个点距离为r的点的总数
[mm,nn]=size(dis);
S2=[];
for r=1:nn
    kk=0;
    s1=0;
    for i1=1:n
        for i2=1:n
            for j1=1:n
                for j2=1:n 
                     minw=min(abs(i1-i2),(8-abs(i1-i2)));
                         minh=min(abs(j1-j2),(8-abs(j1-j2)));                      
                    rr=sqrt(minw^2+minh^2);
                    if (IA1(i1,j1)==1)&&(IA1(i2,j2)==1)&&(rr==dis(r))%统计IA1上所有距离为r的两个值为1的点,加上周期性边界条件。
                     s1=s1+1;
                    else s1=s1+0;
                    end
                end
            end
        end
    end
    if s1~=0
        S2(r)=s1/(n*n*wr(r));%得到相关函数。
    else S2(r)=0;
    end  
end    
end
