function [T,Q,DQ,D2Q]=KorbWodz()
% [T,Q,DQ,D2Q]=KorbWodz
% Rozwi¹zanie zadañ kinematyki o po³o¿eniu, prêdkoœci i przyspieszeniu 
%   dla mechanizmu korbowo-wodzikowego.
% Wyjœcie:
%   T   - tablica do zapisu kolejnych chwil.
%   Q   - tablica do zapisu rozwi¹zañ zadania o po³o¿eniu w kolejnych chwilach.
%   DQ  - tablica do zapisu rozwi¹zañ zadania o prêdkoœci w kolejnych chwilach.
%   D2Q - tablica do zapisu rozwi¹zañ zad. o przyspieszeniu w kolejnych chwilach.
%
%*************************************************************
%   Program stanowi za³¹cznik do ksi¹¿ki:
%   Fr¹czek J., Wojtyra M.: Kinematyka uk³adów wielocz³onowych.
%                           Metody obliczeniowe. WNT 2007.
%   Wersja 1.0
%*************************************************************

% Przybli¿enie startowe (gdy brak rozwi¹zania z poprzedniej chwili)
q=zeros(9,1); dq=zeros(9,1); d2q=zeros(9,1);
lroz=0; % Licznik rozwi¹zañ (s³u¿y do numerowania kolumn w tablicach z wynikami)
dt=0.05; % Odstêp pomiêdzy kolejnymi chwilami

% Rozwi¹zywanie zadañ kinematyki w kolejnych chwilach t
for t=0:dt:1.5
    % Zadanie o po³o¿eniu. 
    % Przybli¿eniem pocz¹tkowym jest rozwi¹zanie z poprzedniej chwili, 
    % powiêkszone o sk³adniki wynikaj¹ce z obliczonej prêdkoœci i przyspieszenia.
    q0=q+dq*dt+0.5*d2q*dt^2;
    q=NewRaph(q0,t); 
    
    dq=Predkosc(q,t);  % Zadanie o prêdkoœci

    d2q=Przyspieszenie(dq,q,t);  % Zadanie o przyspieszeniu

    % Zapis do tablic gromadz¹cych wyniki 
    lroz=lroz+1;
    T(1,lroz)=t; 
    Q(:,lroz)=q;
    DQ(:,lroz)=dq;
    D2Q(:,lroz)=d2q;
end
