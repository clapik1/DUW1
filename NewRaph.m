function q = NewRaph(q, os, ps, wo, wp, t)
%
%   Rozwiazanie ukladu rownan metoda Newtona-Raphsona.
% Wejscie:
%   q0 - przyblizenie startowe,
%   t - chwila, dla ktorej poszukiwane jest rozwiazanie.
% Wyjscie:
%   q - uzyskane rozwiazanie.
%
% Uklad rownan musi byc wpisany w pliku Wiezy.m.
% Macierz Jacobiego ukladu rownan musi byc wpisana w pliku Jakobian.m.
%

F = ones(length(q), 1);
iter = 1;
while norm(F) > 1e-10 && iter < 25
    F = Wiezy(q, os, ps, wo, wp, t);
    Fq = Jakobian(q, os, ps, wo, wp);
    disp(Fq);
    q = q - Fq \ F;
    iter = iter + 1;
end
if iter >= 25
    disp('BLAD: Po 25 iteracjach Newtona-Raphsona nie uzyskano zbieznosci');
end
