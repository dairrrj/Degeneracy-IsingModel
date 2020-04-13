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
        0 0 0 0 0 0 0 0];%��ӦS2t��غ����Ľṹ��
dis=[]; %�������󣬰����п��ܵ����������������
n2 = n/2+1; %�������Եľ���
for ii = 1:n2
   for jj = 1:n2
        rr = sqrt((ii-1)^2 + (jj-1)^2);
           if any(dis == rr) == 0
              dis = [dis,rr];
           end
   end
end
dis = sort(dis);%��С�������С�
[mm,nn] = size(dis);
wr=[];%����һ���ص����������Ϊr�����ص������
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
S2t = s2hs(dis,IA1,wr,n);%��target��غ�����
E=0;
f=exp(1);%����ϵ����
re=[E 0 0];%����ͳ������E���򲢶�g(E)�Ķ��������ִ���H(E)�ľ���
iter=10;
for it = 1:iter %����������
    N=0;
    HEW=0;
    while (HEW==0)
        E1=E;%E1��ʾ��ʼ״̬Ei��
        E2=0;%E2��ʾ��ת״̬Ej��
        ep1 = find( re(:,1)==E1);%����E1�Ƿ�������е�E����ͬ�ģ�����ͬ��ep1������Ӧ��������
        g1 = re(ep1,2);%��ʼ״̬�򲢶ȼ�Ϊg1��

        [x1,y1] = find( IA1(:,:)==1);%�ҵ�����Ϊ1�ĵ��λ��
        [x0,y0] = find( IA1(:,:)==0);%�ҵ�����Ϊ1�ĵ��λ��
      
        l1=length(x1);
        l0=length(x0);

        c1 = ceil( rand*l1);%���ѡȡΪ1�ĵ�
        c0 = ceil( rand*l0);%���ѡȡΪ0�ĵ�

        x11 = x1(c1);
        y11 = y1(c1);

        x00 = x0(c0);
        y00 = y0(c0);

        if (N == 0)&&(it==1)
            S2 = S2t;
        end
        S2tmp = S2;

        S2=s2hss(S2,IA1,dis,wr,n,x11,x00,y11,y00);

        p = IA1(x11,y11);%����1���0��λ�ã����¾���
        IA1(x11,y11) = IA1(x00,y00);
        IA1(x00,y00) = p;
        for r = 1:nn %���¾��������E2��
            E2 = E2+(S2t(r)-S2(r))^2;
        end
        epos = find( re(:,1) == E2);%�ж�E2�Ƿ���ֹ��������֣�epos=1��
        if length(epos) == 0 %E2Ϊ�µ���������û�г��ֹ��ģ�
            E=E2;%��¼E2�������´�ѭ������E2��Ϊ��ʼ״̬��
        re=[re
             E2 0 1];%���¾��󣬽���һ�γ��ֵ�������¼��
            g2=0;%��ת״̬�ļ򲢶ȼ�Ϊg2����E2��һ�γ��֣���Ϊ0��
        else
            g2 = re(epos,2);%��E2�����ֹ���ȡ����¼��gֵ��
        end 

        gg=rand;%�������

        if gg < min(exp(g1-g2),1) && length(epos) ~= 0  %�������Ѿ����ֵ������ҷ��Ϸ�ת���ʣ�
            E=E2;%�´�ѭ������E2��Ϊ��ʼ״̬��
            re(epos,2) = re(epos,2)+log(f);%E2���ڵ�̬�򲢶�����Ϊg(E)=f*g(E)��
            re(epos,3) = re(epos,3)+1;%E2���ֵĴ���+1��
        elseif length(epos) ~= 0 %�������Ѿ����ֵ������������Ϸ�ת���ʣ�
            q = IA1(x11,y11); %���󲻷�ת
            IA1(x11,y11) = IA1(x00,y00);
            IA1(x00,y00) = q;   
            E=E1; %�´�ѭ������E1��Ϊ��ʼ״̬��
            re(ep1,2) = re(ep1,2)+log(f); %E1���ڵ�̬�򲢶�����Ϊg(E)=f*g(E)��
            re(ep1,3) = re(ep1,3)+1; %E1���ֵĴ���+1��
            S2 = S2tmp;
        end

        N=N+1;
        if N>800                %��ֱ��ͼH(E)�Ƿ�ƽ̹��׼
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
    f=sqrt(f);%��������ϵ��f��ʹ��������1��
    if it<iter
        re(:,3)=0;
    end
end
re=sortrows(re);%��re����һ�д�С����
[aa,bb]=size(re);
%��̬�ܶȽ��й�һ��
omegatot =64*63*62*61/(1*2*3*4);%�������ܵ�̬�ܶ�����
lnomegatot = log(omegatot);
minre=min(re(:,2));%��ֹ�������
GTOT=log(sum(exp(re(:,2)-minre)))+minre;%ģ������ܵ�̬�ܶ�����
GN=GTOT-lnomegatot;
re(:,2)=re(:,2)-GN;%���й�һ��
%��̬�ܶ�ͼ
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

subplot(1,2,2)%��ֱ��ͼ
bar(re(:,3))

%����ͳ������E�ͼ򲢶�g(E),(��һ��״̬jj+1�Ļ���)��

toc
