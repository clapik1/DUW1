function d2q=Przyspieszenie(dq,q,t)
% d2q=Przyspieszenie(dq,q,t)
%   Procedura do rozwi¹zywania zadania o przyspieszeniu.
%   Zadania o po³o¿eniu i prêdkoœci musz¹ byæ rozwi¹zane wczeœniej.
% Wejœcie:
%   dq - pochodne wspó³rzêdnych absolutnych wzglêdem czasu,
%   q  - wspó³rzêdne absolutne uk³adu wielocz³onowego,
%   t  - chwila, dla której poszukiwane jest rozwi¹zanie.
% Wyjœcie:
%   d2q - obliczone drugie pochodne wspó³rzêdnych absolutnych wzglêdem czasu.
%
%*************************************************************
%   Program stanowi za³¹cznik do ksi¹¿ki:
%   Fr¹czek J., Wojtyra M.: Kinematyka uk³adów wielocz³onowych.
%                           Metody obliczeniowe. WNT 2007.
%   Wersja 1.0
%*************************************************************

Dane; % Wczytanie danych o wymiarach mechanizmu i usytuowaniu par kinematycznych

% Przypisanie elementom wektora q czytelnych nazw
r1=q(1:2); fi1=q(3);   r2=q(4:5); fi2=q(6);  r3=q(7:8); fi3=q(9);

% Obliczenie macierzy kosinusów kierunkowych
Rot1=Rot(fi1);  Rot2=Rot(fi2);  Rot3=Rot(fi3);

% Przypisanie elementom wektora dq czytelnych nazw
dfi1=dq(3);   dfi2=dq(6);   dfi3=dq(9);

% Obliczenie wektora prawych stron
gam(1:2,1)=-Rot1*sB1*dfi1^2;
gam(3:4,1)=Rot1*sA1*dfi1^2-Rot2*sB2*dfi2^2;
gam(5:6,1)=Rot2*sA2*dfi2^2-Rot3*sB2*dfi2^2;
gam(7,1)=0;
gam(8,1)=-v0'*Rot3*sA3*dfi3^2;
gam(9,1)=2;

% Obliczenie macierzy uk³adu równañ
Fq=Jakobian(q);

% Obliczenie przyspieszenia (rozwi¹zanie uk³adu równañ)
d2q=Fq\gam; 
