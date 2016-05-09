function F=Wiezy(q,t)
% F=Wiezy(q,t)
%   Procedura wspolpracujaca z NewRaph.
%   Sluzy do obliczania wartosci funkcji opisujacych wiezy.
% Wejscie:
%   q - wspolrzedne absolutne ukladu wieloczlonowego,
%   t - aktualna chwila.
% Wyjscie:
%   F - wartosci funkcji.
%

F = zeros(length(q0),1);

rozos = size(os);
rozF = 1;
for i=1:rozos(1) % petla po wszystkich parach obrotowych
    [ri,~,Roti] = FromQ(q,os(i,1));
    [rj,~,Rotj] = FromQ(q,os(i,2));
    sA = os(i,3:4)';
    sB = os(i,5:6)';
    F(rozF:rozF+1,1) = ri + Roti * sA - (rj + Rotj * sB); % wzor 2.16 na pare obrotowa
    rozF = rozF + 2;
end

rozps = size(ps);
for i=1:rozps(1) % petla po wszystkich parach postepowych
    [ri,fii,Roti] = FromQ(q,ps(i,1));
    [rj,fij,Rotj] = FromQ(q,ps(i,2));
    fi0 = ps(i,3);
    v = ps(i,4:5)';
    sA = ps(i,6:7)';
    sB = ps(i,8:9)';
    F(rozF,1) = fii - fij - fi0; % wzor 2.17 na warunek kata w parze postepowej
    F(rozF+1,1) = (Rotj * v)'*(rj - ri - Roti * sA) + v' * sB; % wzor 2.20 na brak ruchu w kierunku v
    rozF = rozF + 2;
end

rozwo = size(wo);
for i=1:rozwo % petla po wszystkich wymuszeniach w parach obrotowych
    [~,fii,~] = FromQ(q,os(rozwo(i,1),1));
    [~,fij,~] = FromQ(q,os(rozwo(i,1),2));
    F(rozF,1) = fii - fij - Wymuszenie(rozwo(i,2),t); % wzor 2.25 na wzajemna orientacje czlonow
    rozF = rozF + 1;
end

rozwp = size(wp);
for i=1:rozwp(1) % petla po wszystkich wymuszeniach w parach postepowych
    [ri,~,Roti] = FromQ(q,ps(rozwp(i,1),1));
    [rj,~,Rotj] = FromQ(q,ps(rozwp(i,1),2));
    u = ps(rozwp(i,1),4:5)';
    sA = ps(rozwp(i,1),6:7)';
    sB = ps(rozwp(i,1),8:9);
    % wzor 2.26 na wymuszenie w parze postepowej
    F(rozF,1) = (Rotj * u)'*(rj + Rotj * sB - ri - Roti * sA) - Wymuszenie(rozwp(i,2),t);
    rozF = rozF + 1;
end
