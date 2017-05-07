function [T, Q, DQ, D2Q, EXPQ] = DUW1()
%
% Rozwiazanie zadan kinematyki o polozeniu, predkosci i przyspieszeniu 
%   dla mechanizmu korbowo-wodzikowego.
% Wyjscie:
%   T   - tablica do zapisu kolejnych chwil.
%   Q   - tablica do zapisu rozwiazan zadania o polozeniu w kolejnych chwilach.
%   DQ  - tablica do zapisu rozwiazan zadania o predkosci w kolejnych chwilach.
%   D2Q - tablica do zapisu rozwiazan zad. o przyspieszeniu w kolejnych chwilach.
%

[q, os, ps, wo, wp, de] = Wczyt();

Om = [0 -1; 1 0];
lroz = 0;
dt = 0.05;
dq = zeros(length(q), 1);
d2q = zeros(length(q), 1);

s = 30 / dt + 1;
T = zeros(1, s);
Q = zeros(length(q), s);
DQ = zeros(length(q), s);
D2Q = zeros(length(q), s);
EXPQ = zeros(size(de, 1) * 6, s);

for t = 0:dt:30
    % Przyblizeniem poczatkowym jest rozwiazanie z poprzedniej chwili, 
    % powiekszone o skladniki wynikajace z obliczonej predkosci i przyspieszenia.
    q = q + dq * dt + 0.5 * d2q * dt^2;
    q = NewRaph(q, os, ps, wo, wp, t); 
    dq = Predkosc(q, os, ps, wo, wp, t);
    d2q = Przyspieszenie(dq, q, os, ps, wo, wp, t);
    
    %{
    hold on;
    fir = zeros(length(q)/3, 2);
    tab = zeros(length(q)/3, 2);
    for i = 1:size(os,1)
        if os(i,1) ~= 0
            w1 = (Rot(q(os(i,1)*3)) * os(i,3:4)')';
            if tab(os(i,1),1) ~= 0 && tab(os(i,1),2) ~= 0
                plot([tab(os(i,1),1); w1(1)], [tab(os(i,1),2); w1(1)]);
            else
                fir(os(i,1),:) = w1;
            end
            tab(os(i,1),:) = w1;
        end
        if os(i,2) ~= 0
            w1 = (Rot(q(os(i,2)*3)) * os(i,5:6)')';
            if tab(os(i,2),1) ~= 0 && tab(os(i,2),2) ~= 0
                plot([tab(os(i,2),1); w1(1)], [tab(os(i,2),2); w1(1)]);
            else
                fir(os(i,2),:) = w1;
            end
            tab(os(i,2),:) = w1;
        end
    end
    for i = 1:size(tab,1)
        plot([tab(1); fir(1)], [tab(2); fir(2)]);
    end
    hold off;
    %}
    
    % Zapis do tablic gromadzacych wyniki 
    lroz = lroz + 1;
    T(:, lroz) = t; 
    Q(:, lroz) = q;
    DQ(:, lroz) = dq;
    D2Q(:, lroz) = ddq;
    
    for i = 1:size(de, 1)
        sA = de(i,2:3)';
        [ri,~,Roti] = FromQ(q, de(i, 1));
        [dri,dfii,~] = FromQ(dq, de(i, 1));
        [ddri,ddfii,~] = FromQ(ddq, de(i, 1));
        EXPQ(i*6-5:i*6-4, lroz) = ri + Roti * sA;
        EXPQ(i*6-3:i*6-2, lroz) = dri + Om * Roti * sA * dfii;
        EXPQ(i*6-1:i*6, lroz) = ddri + Om * Roti * sA * ddfii - Roti * sA * dfii^2;
    end
end
