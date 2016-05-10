function [T, Q, DQ, D2Q] = DUW1()
%
% Rozwiazanie zadan kinematyki o polozeniu, predkosci i przyspieszeniu 
%   dla mechanizmu korbowo-wodzikowego.
% Wyjscie:
%   T   - tablica do zapisu kolejnych chwil.
%   Q   - tablica do zapisu rozwiazan zadania o polozeniu w kolejnych chwilach.
%   DQ  - tablica do zapisu rozwiazan zadania o predkosci w kolejnych chwilach.
%   D2Q - tablica do zapisu rozwiazan zad. o przyspieszeniu w kolejnych chwilach.
%

[q, os, ps, wo, wp] = Wczyt();

lroz = 0;
dt = 0.05;
dq = zeros(length(q), 1);
d2q = zeros(length(q), 1);

s = 1.5 / dt + 1;
T = zeros(1, s);
Q = zeros(length(q), s);
DQ = zeros(length(q), s);
D2Q = zeros(length(q), s);

for t = 0:dt:1.5
    % Przyblizeniem poczatkowym jest rozwiazanie z poprzedniej chwili, 
    % powiekszone o skladniki wynikajace z obliczonej predkosci i przyspieszenia.
    q = q + dq * dt + 0.5 * d2q * dt^2;
    q = NewRaph(q, os, ps, wo, wp, t); 
    dq = Predkosc(q, os, ps, wo, wp, t);
    d2q = Przyspieszenie(dq, q, os, ps, wo, wp, t);

    % Zapis do tablic gromadzacych wyniki 
    lroz = lroz + 1;
    T(1, lroz) = t; 
    Q(:, lroz) = q;
    DQ(:, lroz) = dq;
    D2Q(:, lroz) = d2q;
end
