function F = Wiezy(q, os, ps)
%
%   Procedura wspolpracujaca z NewRaph.
%   Sluzy do obliczania wartosci funkcji opisujacych wiezy.
% Wejscie:
%   q - wspolrzedne absolutne ukladu wieloczlonowego,
%   t - aktualna chwila.
% Wyjscie:
%   F - wartosci funkcji.
%

F = zeros(2*(size(os,1)+size(ps,1)),1);
rozF = 1;

for i=1:size(os, 1) % petla po wszystkich parach obrotowych
    [ri,~,Roti] = FromQ(q,os(i,1));
    [rj,~,Rotj] = FromQ(q,os(i,2));
    sA = os(i,3:4)';
    sB = os(i,5:6)';
    F(rozF:rozF+1,1) = ri + Roti * sA - (rj + Rotj * sB); % wzor 2.16 na pare obrotowa
    rozF = rozF + 2;
end

for i=1:size(ps, 1) % petla po wszystkich parach postepowych
    [ri,fii,Roti] = FromQ(q,ps(i,1));
    [rj,fij,Rotj] = FromQ(q,ps(i,2));
    fi0 = ps(i,3);
    u = ps(i,4:5)';
    v = [-u(2); u(1)];
    sA = ps(i,6:7)';
    sB = ps(i,8:9)';
    F(rozF,1) = fii - fij - fi0; % wzor 2.17 na warunek kata w parze postepowej
    F(rozF+1,1) = (Rotj * v)'*(rj - ri - Roti * sA) + v' * sB; % wzor 2.20 na brak ruchu w kierunku v
    rozF = rozF + 2;
end

% for i=1:size(wo, 1) % petla po wszystkich wymuszeniach w parach obrotowych
%     [~,fii,~] = FromQ(q,os(wo(i,1),1));
%     [~,fij,~] = FromQ(q,os(wo(i,1),2));
%     F(rozF,1) = fii - fij - Wymuszenie(wo(i,2),t); % wzor 2.25 na wzajemna orientacje czlonow
%     rozF = rozF + 1;
% end
% 
% for i=1:size(wp, 1) % petla po wszystkich wymuszeniach w parach postepowych
%     [ri,~,Roti] = FromQ(q,ps(wp(i,1),1));
%     [rj,~,Rotj] = FromQ(q,ps(wp(i,1),2));
%     u = ps(wp(i,1),4:5)';
%     sA = ps(wp(i,1),6:7)';
%     sB = ps(wp(i,1),8:9)';
%     % wzor 2.26 na wymuszenie w parze postepowej
%     F(rozF,1) = (Rotj * u)'*(rj + Rotj * sB - ri - Roti * sA) - Wymuszenie(wp(i,2),t);
%     rozF = rozF + 1;
% end
