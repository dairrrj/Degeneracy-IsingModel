function S2=s2hs(dis,IA1,wr,n)
%disΪ����r����
%IA1ΪҪ����غ����Ľṹ�ı�ʾ����
%wr��ʾ����һ�������Ϊr�ĵ������
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
                    if (IA1(i1,j1)==1)&&(IA1(i2,j2)==1)&&(rr==dis(r))%ͳ��IA1�����о���Ϊr������ֵΪ1�ĵ�,���������Ա߽�������
                     s1=s1+1;
                    else s1=s1+0;
                    end
                end
            end
        end
    end
    if s1~=0
        S2(r)=s1/(n*n*wr(r));%�õ���غ�����
    else S2(r)=0;
    end  
end    
end
