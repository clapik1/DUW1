function R=Rot(fi)
% R=Rot(fi)
%   Obliczanie macierzy rotacji - przypadek p³aski (2x2)
% Wejœcie:
%   fi - kat fi (w radianach).
% Wyjœcie:
%   R - macierz kosinusów kierunkowych.
%
%*************************************************************
%   Program stanowi za³¹cznik do ksi¹¿ki:
%   Fr¹czek J., Wojtyra M.: Kinematyka uk³adów wielocz³onowych.
%                           Metody obliczeniowe. WNT 2007.
%   Wersja 1.0
%*************************************************************

R=[cos(fi) -sin(fi);sin(fi) cos(fi)];

