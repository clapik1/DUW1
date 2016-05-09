function dq=Predkosc(q,t)
% dq=Predkosc(q,t)
%   Procedura do rozwi¹zywania zadania o prêdkoœci.
%   Zadanie o po³o¿eniu musi byæ rozwi¹zane wczeœniej.
% Wejœcie:
%   q - wspó³rzêdne absolutne uk³adu wielocz³onowego,
%   t - chwila, dla której poszukiwane jest rozwi¹zanie.
% Wyjœcie:
%   dq - obliczone pochodne wspó³rzêdnych absolutnych wzglêdem czasu.
%
%*************************************************************
%   Program stanowi za³¹cznik do ksi¹¿ki:
%   Fr¹czek J., Wojtyra M.: Kinematyka uk³adów wielocz³onowych.
%                           Metody obliczeniowe. WNT 2007.
%   Wersja 1.0
%*************************************************************

% Prawa strona uk³adu równañ liniowych - pochodne wiêzów wzglêdem czasu.
Ft=[zeros(8,1);-2*t];

% Macierz uk³adu równañ
Fq=Jakobian(q);

% Obliczenie prêdkoœci
dq=-Fq\Ft;
