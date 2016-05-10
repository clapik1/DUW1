function [T, Q, DQ, D2Q] = DUW1_korbwodz()
%
% Rozwiazanie zadan kinematyki o polozeniu, predkosci i przyspieszeniu 
%   dla mechanizmu korbowo-wodzikowego.
% Wyjscie:
%   T   - tablica do zapisu kolejnych chwil.
%   Q   - tablica do zapisu rozwiazan zadania o polozeniu w kolejnych chwilach.
%   DQ  - tablica do zapisu rozwiazan zadania o predkosci w kolejnych chwilach.
%   D2Q - tablica do zapisu rozwiazan zad. o przyspieszeniu w kolejnych chwilach.
%

% [q, os, ps, wo, wp] = Wczyt();
a=1; b=5; c=4;  % Wymiary mechanizmu
q = [0; 0; pi/2; 0; 1; atan(3/4); 4; 4; 0];
os = [0 1 0 0 0 0; 1 2 a 0 0 0; 2 3 b 0 0 0];
ps = [0 3 0 1 0 0 c 0 0];
%wo = [1 4];
wo = [];
%wp = [];
wp = [1 5];

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
    T(:, lroz) = t; 
    Q(:, lroz) = q;
    DQ(:, lroz) = dq;
    D2Q(:, lroz) = d2q;

    w1 = Rot(q(3)) * os(1,5:6)';
    w2 = Rot(q(3)) * os(2,3:4)';
    plot([w1(1) w2(1)] + q(1), [w1(2), w2(2)] + q(2), 'r');
    hold on
    w1 = Rot(q(6)) * os(2,5:6)';
    w2 = Rot(q(6)) * os(3,3:4)';
    plot([w1(1) w2(1)] + q(4), [w1(2), w2(2)] + q(5), 'b');
    hold off
    
    axis([-5 5 -5 5]);
    pause(0.1);
end
