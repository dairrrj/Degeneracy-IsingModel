clear;
clc;
tic
n=8;
IA1=[0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0
        0 0 1 0 0 1 0 0
        0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 
        0 0 1 0 0 1 0 0
        0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0];%对应S2t相关函数的结构。
dis=[]; %求距离矩阵，把所有可能的两点间距离存下来。
n2 = n/2+1; %有周期性的距离
for ii = 1:n2
   for jj = 1:n2
        rr = sqrt((ii-1)^2 + (jj-1)^2);
           if any(dis == rr) == 0
              dis = [dis,rr];
           end
   end
end
dis = sort(dis);%由小到大排列。
[mm,nn] = size(dis);
wr=[];%任意一像素点距离它长度为r的像素点的总数
for ww = 1:nn
    kk=0;
    for ii = 1:n
        for jj = 1:n
            minw = min(abs(ii-1) , (n-abs(ii-1)));
            minh = min(abs(jj-1) , (n-abs(jj-1)));
            rr = sqrt(minw^2+minh^2);           
            if rr == dis(ww)
                kk = kk+1;
            end   
        end
    end
    wr(ww) = kk;
end
S2t = s2hs(dis,IA1,wr,n);%求target相关函数。
E=0;
f=exp(1);%修正系数。
re=[E 0 0];%生成统计能量E、简并度g(E)的对数、出现次数H(E)的矩阵。
iter=10;
for it = 1:iter %迭代次数；
    N=0;
    HEW=0;
    while (HEW==0)
        E1=E;%E1表示初始状态Ei；
        E2=0;%E2表示翻转状态Ej。
        ep1 = find( re(:,1)==E1);%查找E1是否与矩阵中的E有相同的，若相同，ep1等于相应的行数；
        g1 = re(ep1,2);%初始状态简并度记为g1。

        [x1,y1] = find( IA1(:,:)==1);%找到所有为1的点的位置
        [x0,y0] = find( IA1(:,:)==0);%找到所有为1的点的位置
      
        l1=length(x1);
        l0=length(x0);

        c1 = ceil( rand*l1);%随机选取为1的点
        c0 = ceil( rand*l0);%随机选取为0的点

        x11 = x1(c1);
        y11 = y1(c1);

        x00 = x0(c0);
        y00 = y0(c0);

        if (N == 0)&&(it==1)
            S2 = S2t;
        end
        S2tmp = S2;

        S2=s2hss(S2,IA1,dis,wr,n,x11,x00,y11,y00);

        p = IA1(x11,y11);%交换1点和0点位置，更新矩阵
        IA1(x11,y11) = IA1(x00,y00);
        IA1(x00,y00) = p;
        for r = 1:nn %求新矩阵的能量E2。
            E2 = E2+(S2t(r)-S2(r))^2;
        end
        epos = find( re(:,1) == E2);%判断E2是否出现过，若出现，epos=1。
        if length(epos) == 0 %E2为新的能量，即没有出现过的；
            E=E2;%记录E2，并在下次循环中由E2作为初始状态；
        re=[re
             E2 0 1];%更新矩阵，将第一次出现的能量记录；
            g2=0;%翻转状态的简并度记为g2，若E2第一次出现，记为0；
        else
            g2 = re(epos,2);%若E2曾出现过，取曾记录的g值。
        end 

        gg=rand;%随机数。

        if gg < min(exp(g1-g2),1) && length(epos) ~= 0  %遍历到已经出现的能量且符合翻转概率；
            E=E2;%下次循环中由E2作为初始状态；
            re(epos,2) = re(epos,2)+log(f);%E2所在的态简并度修正为g(E)=f*g(E)；
            re(epos,3) = re(epos,3)+1;%E2出现的次数+1。
        elseif length(epos) ~= 0 %遍历到已经出现的能量但不符合翻转概率；
            q = IA1(x11,y11); %矩阵不翻转
            IA1(x11,y11) = IA1(x00,y00);
            IA1(x00,y00) = q;   
            E=E1; %下次循环中由E1作为初始状态；
            re(ep1,2) = re(ep1,2)+log(f); %E1所在的态简并度修正为g(E)=f*g(E)；
            re(ep1,3) = re(ep1,3)+1; %E1出现的次数+1。
            S2 = S2tmp;
        end

        N=N+1;
        if N>800                %求直方图H(E)是否平坦标准
            [aa,bb]=size(re);
            ttt=0;
            for ee=1:aa
                ttt=ttt+re(ee,3);
            end
            tttt=ttt/aa*(0.80);
            if min(re(:,3))>=tttt
                HEW=1;
                N
            else
                HEW=0;
            end
        end
    end
    f=sqrt(f);%修正修正系数f，使其逐步趋于1；
    if it<iter
        re(:,3)=0;
    end
end
re=sortrows(re);%将re按第一列大小排列
[aa,bb]=size(re);
%对态密度进行归一化
omegatot =64*63*62*61/(1*2*3*4);%理论上总的态密度数量
lnomegatot = log(omegatot);
minre=min(re(:,2));%防止溢出处理
GTOT=log(sum(exp(re(:,2)-minre)))+minre;%模拟出的总的态密度数量
GN=GTOT-lnomegatot;
re(:,2)=re(:,2)-GN;%进行归一化
%作态密度图
NGEE=re(:,2);
NGE = zeros(aa,1);
for e=1:aa
    for kkk=1:aa
        if re(kkk,1)<=re(e,1)
            NGE(e)=NGE(e)+exp(NGEE(kkk));
        end
    end
end
sizere=size(re);
ladder=[];
for i=1:sizere(1)-1
    ladder=[ladder
            re(i,1) NGE(i)
            re(i+1,1) NGE(i)];
end
ladder=[ladder
        re(sizere(1),1) NGE(sizere(1))];
ladder(1,1) = 1e-6;

subplot(1,2,1)
loglog(ladder(:,1),ladder(:,2),'.-')

subplot(1,2,2)%作直方图
bar(re(:,3))

%返回统计能量E和简并度g(E),(下一个状态jj+1的基础)；

toc
