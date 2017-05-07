function [t, q, dq, ddq]=Calkuj_Baumgarte(q0,tK,os,ps,ws,wst,m)
% function [t, Y]=Calkuj_Baumgarte(q0,tK,os,ps,ws,wst,m)

%Ca³kowanie rónañ ruchu metod¹ Rungego-Kutty - ode45
%Stabilizacja wiêzów metod¹ Baumgarte'a

tic %Pocz¹tek pomiaru czasu wykonywania obliczeñ
%Przedzia³ ca³kowania
t0=0; %t0 - pocz¹tek tK - koniec

%Warunki pocz¹tkowe
% q0=[0 0 0 1/sqrt(2) -1/sqrt(2) pi/4]';
dq0=zeros(size(q0,1),1);
Y0=[q0;dq0];

%Ca³kowanie
[t,Y]=ode45(@(t,Y) H(t,Y,os,ps,ws,wst,m),[t0 tK],Y0);

ddq=zeros(size(Y,1), size(Y,2));
for iter=1:size(Y,1)-1
   ddq(iter, :)=H(t(iter), Y(iter,:)', os, ps, ws, wst, m); 
end

ddq=ddq(:, size(q0,1)+1:end)';
q=Y(:, 1:size(q0,1))';
dq=Y(:, size(q0,1)+1:end)';

toc %Koniec pomiary czasu wykonywania obliczeñ

%Drukowanie wyników dla czasu równego tK
n=size(Y,1);
qK=Y(n,1:size(q0,1))'; dqK=Y(n,(size(q0,1)+1):(2*size(q0,1)))';     %Po³o¿enia i prêdkoœci w chwili tK
q__dq=[qK dqK];                                                      %Drukowanie
[F]=Macierze(qK,dqK,os,ps,ws,wst,m);                                %Obliczanie lewej strony równañ wiêzów
norm_F=norm(F)                                                      %Miara rozerwania wiêzów po czasie tK

end

function [dY]= H(t,Y,os,ps,ws,wst,m)
%Funkcja podca³kowa
alf=5; bet=5;                           % Wspó³czynniki odpowiadaj¹ce za stabilizacjê met. Baumgarte'a
q=Y(1:(size(Y,1)/2), :);                % czytelne nazwy po³o¿eñ i prêdkoœci
dq=Y(((size(Y,1)/2)+1):size(Y,1),:);
[F,Fq,G,M,Q]=Macierze(q,dq,os,ps,ws,wst,m);            
A=[M,Fq' ;Fq,zeros(size(Fq,1),size(Fq,1))];               
b=[Q;G-2*alf*Fq*dq-bet^2*F];           

x=A\b;                                  % Obliczenie przyspieszeñ i mno¿ników Lagrange'a
dY(1:(size(Y,1)/2),1)=dq;               % Obliczenie wektora prawych stron równania ró¿niczkowego

dY(((size(Y,1)/2)+1):size(Y,1),1)=x(1:(size(Y,1)/2),1);
end
