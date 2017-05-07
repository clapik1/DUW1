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