function Fq=Jakobian(q)
% Fq=Jakobian(q)
%   Procedura wspó³pracuj¹ca z NewRaph.
%   S³u¿y do obliczania macierzy Jacobiego równañ wiêzów.
% Wejœcie:
%   q - wspó³rzêdne absolutne uk³adu wielocz³onowego.
% Wyjœcie:
%   Fq - obliczona macierz Jacobiego.
%
%*************************************************************
%   Program stanowi za³¹cznik do ksi¹¿ki:
%   Fr¹czek J., Wojtyra M.: Kinematyka uk³adów wielocz³onowych.
%                           Metody obliczeniowe. WNT 2007.
%   Wersja 1.0
%*************************************************************

Dane; % Wczytanie danych o wymiarach mechanizmu i usytuowaniu par kinematycznych

Om=[0 -1;1 0];  % Sta³a macierz

% Przypisanie elementom wektora q czytelnych nazw
r1=q(1:2); fi1=q(3);   r2=q(4:5); fi2=q(6);  r3=q(7:8); fi3=q(9);

% Obliczenie macierzy kosinusów kierunkowych
Rot1=Rot(fi1);  Rot2=Rot(fi2);  Rot3=Rot(fi3);

% Macierz Jaocobiego - pocz¹tkowo zerowa
Fq=zeros(9,9);
% Macierz Jaocobiego - wype³nianie niezerowych elementów
Fq(1:2,1:2)=-eye(2);
Fq(1:2,3)=-Om*Rot1*sB1;

Fq(3:4,1:2)=eye(2);
Fq(3:4,3)=Om*Rot1*sA1;
Fq(3:4,4:5)=-eye(2);
Fq(3:4,6)=-Om*Rot2*sB2;

Fq(5:6,4:5)=eye(2);
Fq(5:6,6)=Om*Rot2*sA2;
Fq(5:6,7:8)=-eye(2);
Fq(5:6,9)=-Om*Rot3*sB3;

Fq(7,9)=1;

Fq(8,7:8)=-v0';
Fq(8,9)=-v0'*Om*Rot3*sA3;

Fq(9,3)=1; 
