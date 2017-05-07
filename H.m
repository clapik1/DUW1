function [dY]= H(t,Y,os,ps,ws,wst,m)
%Funkcja podca�kowa
alf=5; bet=5;                           % Wsp�czynniki odpowiadaj�ce za stabilizacj� met. Baumgarte'a
q=Y(1:(size(Y,1)/2), :);                % czytelne nazwy po�o�e� i pr�dko�ci
dq=Y(((size(Y,1)/2)+1):size(Y,1),:);
[F,Fq,G,M,Q]=Macierze(q,dq,os,ps,ws,wst,m);            
A=[M,Fq' ;Fq,zeros(size(Fq,1),size(Fq,1))];               
b=[Q;G-2*alf*Fq*dq-bet^2*F];           

x=A\b;                                  % Obliczenie przyspiesze� i mno�nik�w Lagrange'a
dY(1:(size(Y,1)/2),1)=dq;               % Obliczenie wektora prawych stron r�wnania r�niczkowego

dY(((size(Y,1)/2)+1):size(Y,1),1)=x(1:(size(Y,1)/2),1);
end