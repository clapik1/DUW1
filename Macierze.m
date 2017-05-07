function [F,Fq,G,M,Q]=Macierze(q,dq,os,ps,ws,wst,m)
% %[F,Fq,G,M,Q]=Macierze(q,dq)
% %
% %Wejœcie:
% % q(6,1)  - wektor wspó³rzêdnych absolutnych
% % dq(6,1) - wektor pochodnych wspó³rzêdnych absolutnych
% %
% %Wyjœcie:
% % F(4,1)  - lewe strony równañ wiêzów
% % Fq(4,6) - macierz Jacobiego
% % G(4,1)  - wektor gamma
% % M(6,6)  - macierz masowa
% % Q(6,1)  - wektor si³ uogólnionych
% 
% 
% % Wymiary mechanizmu:
% a=1;                                     % D³ugoœæ wahad³a
% f0=0; v0=[0;1]; sA1=[0;0]; sB0=[100;0];  % Para postêpowa 1-0
% sB2=[0;a];                               % Para obrotowa 1-2 (sA1 zdefiniowane wczeœniej)
% rE0=[4*a;-a]; sD2=[0;-a];                % Punkty mocowania t³umika
% % Sta³a macierz:
% Om=[0 -1;1 0]; 
% Przyspieszenie ziemskie (wektor):
g=[0; -9.81];
Om = [0 -1; 1 0];
% % Wspó³czynnik t³umienia:
% c=10;
% % Masy mechanizmu:
% m1=5; J1=2; m2=3; J2=m2*a^2/3;
% 
% % Przypisanie elementom wektora q czytelnych nazw
% r1=q(1:2); fi1=q(3);   r2=q(4:5); fi2=q(6);
% % Przypisanie elementom wektora dq czytelnych nazw
% dr1=dq(1:2,1);  dfi1=dq(3,1);   dr2=dq(4:5,1);  dfi2=dq(6,1);
% 
% % Obliczenie macierzy kosinusów kierunkowych
% Rot1=Rot(fi1);  Rot2=Rot(fi2); 
% 
% % Lewe strony równañ wiêzów
% F(1,1)=fi1-f0;
% F(2,1)=v0'*(sB0-r1-Rot1*sA1);
% F(3:4,1)=r1+Rot1*sA1-(r2+Rot2*sB2);
% 
% % Macierz Jacobiego
% Fq=zeros(4,6);
% Fq(1,3)=1;
% Fq(2,1:2)=-v0';
% Fq(2,3)=-v0'*Om*Rot1*sA1;
% Fq(3:4,1:2)=eye(2);
% Fq(3:4,3)=Om*Rot1*sA1;
% Fq(3:4,4:5)=-eye(2);
% Fq(3:4,6)=-Om*Rot2*sB2;
% 
% % Wektor Gamma
% G(1,1)=0;
% G(2,1)=-v0'*Rot1*sA1*dfi1^2;
% G(3:4,1)=Rot1*sA1*dfi1^2-Rot2*sB2*dfi2^2;

% % Lewe strony równañ wiêzów
F = Wiezy(q,os,ps);

% % Macierz Jacobiego
Fq = Jakobian(q, os, ps);

% % Wektor Gamma
G = Gamma(dq, q, os, ps);

% Macierz masowa
M=zeros(size(m, 1)*3, size(m, 1)*3);
for i=1:size(m, 1)
    M(3*i-2:3*i, 3*i-2:3*i)=diag([m(i, 1), m(i, 1), m(i, 2)]);
end

% Sily uogolnione - grawitacja
Q=zeros(size(m, 1)*3, 1);
for i=1:size(m, 1)
    Q(3*i-2:3*i, 1)=[m(i, 1)*g; 0];
end

% Sily uogolnione - przylozone
for k=1:size(ws, 1)
    i=ws(k,1);
    s=ws(k, 2:3)';
    [~, ~, Roti] = FromQ(q, i);
    F=ws(k, 4:5)';
    
    Q(3*i-2:3*i, 1) = Q(3*i-2:3*i, 1) + [eye(2); (Om*Roti*s)']*F;
end

% Si³y uogólnione - od elementów sprê¿ysto t³umi¹cych
for k=1:size(wst, 1)
 
    %Element t³umi¹cy
    i=wst(k,1);
    j=wst(k,2);
    sA=[wst(k,3),wst(k,4)]';                                            %wspó³rzêdne wktora punktu przy³o¿enia si³a i
    sB=[wst(k,5),wst(k,6)]';                                            %wspó³rzêdne wktora punktu przy³o¿enia si³a j
    
    
    [ri, ~, Roti] = FromQ(q, i);
    [rj, ~, Rotj] = FromQ(q, j);
    d= rj + Rotj*sB - ri - Roti*sA;                                     %wzór b str. 10
    u=d/norm(d);                                                        %wzór c

    [dri, dfii, ~] = FromQ(dq, i);
    [drj, dfij, ~] = FromQ(dq, j);
    
    dd=u'*(drj + Om*Rotj*sB*dfij - dri - Om*Roti*sA*dfii);              %wzór d
    
    F=u*wst(k,9)*dd;                                                    %wzór e
    
    if(i~=0)
        Q(3*i-2:3*i, 1) = Q(3*i-2:3*i, 1) + [eye(2); (Om*Roti*sA)']*F;      %wzór f
    end
    
    if(j~=0)
        Q(3*j-2:3*j, 1) = Q(3*j-2:3*j, 1) + [eye(2); (Om*Rotj*sB)']*(-F);
    end
    
    %Element sprê¿ysty
    l0=wst(k,8);                                                        %d³ugoœæ swobodna si³ownika
    F=u*wst(k,7)*(norm(d)-l0);                                            %si³a sprê¿ystoœci (minus uwzglêdnony w przeciwnym zwrocie wersora u)
    
    if(i~=0)
        Q(3*i-2:3*i, 1) = Q(3*i-2:3*i, 1) + [eye(2); (Om*Roti*sA)']*F;    
    end
    
    if(j~=0)
        Q(3*j-2:3*j, 1) = Q(3*j-2:3*j, 1) + [eye(2); (Om*Rotj*sB)']*(-F);    
        %( trzecim wierszu wstawi³em dwa zera, bo sprê¿yny s¹ zawsze przy³¹czone do œrodkó mas, chyba?)
    end
end

% % Lewe strony równañ wiêzów
F = Wiezy(q,os,ps);

% Si³y uogólnione - przylo¿one
% d=rE0-r2-Rot2*sD2;
% u=d/sqrt(d'*d);
% dd=-u'*(dr2+Om*Rot2*sD2*dfi2);
% Q(4:6,1)=Q(4:6,1)+[eye(2);(Om*Rot2*sD2)']*u*c*dd;