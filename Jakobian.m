function Fq = Jakobian(q, os, ps, wo, wp)
%
%   Procedura wspolpracujaca z NewRaph.
%   Sluzy do obliczania macierzy Jacobiego rownan wiezow.
% Wejscie:
%   q - wspolrzedne absolutne ukladu wieloczlonowego.
% Wyjscie:
%   Fq - obliczona macierz Jacobiego.
%

Om = [0 -1; 1 0];
Fq = zeros(length(q),length(q));
rozF = 1;

for i=1:size(os, 1) % petla po wszystkich parach obrotowych
    sA = os(i,3:4)';
    sB = os(i,5:6)';
    [~,~,Roti] = FromQ(q,os(i,1));
    [~,~,Rotj] = FromQ(q,os(i,2));
    if os(i,1) ~= 0 % wazne tylko, gdy czlon nie jest podstawa
        Fq(rozF:rozF+1,os(i,1)*3-2:os(i,1)*3-1) = eye(2); % 2.29
        Fq(rozF:rozF+1,os(i,1)*3) = Om * Roti * sA; % 2.30
    end
    if os(i,2) ~= 0
        Fq(rozF:rozF+1,os(i,2)*3-2:os(i,2)*3-1) = -eye(2); % 2.31
        Fq(rozF:rozF+1,os(i,2)*3) = -Om * Rotj * sB; % 2.32
    end
    rozF = rozF + 2;
end

for i=1:size(ps, 1) % petla po wszystkich parach postepowych
    u = ps(i,4:5)';
    v = [-u(2); u(1)];
    sA = ps(i,6:7)';
    [ri,~,Roti] = FromQ(q,ps(i,1));
    [rj,~,Rotj] = FromQ(q,ps(i,2));
    if ps(i,1) ~= 0
        Fq(rozF,ps(i,1)*3) = 1; % 2.42
        Fq(rozF+1,ps(i,1)*3-2:ps(i,1)*3-1) = -(Rotj * v)'; % 2.47
        Fq(rozF+1,ps(i,1)*3) = -(Rotj * v)' * Om * Roti * sA; % 2.48
    end
    if ps(i,2) ~= 0
        Fq(rozF,ps(i,2)*3) = -1; % 2.44
        Fq(rozF+1,ps(i,2)*3-2:ps(i,2)*3-1) = (Rotj * v)'; % 2.49
        Fq(rozF+1,ps(i,2)*3) = -(Rotj * v)' * Om * (rj - ri - Roti * sA); % 2.50
    end
    rozF = rozF + 2;
end

for i=1:size(wo, 1) % petla po wszystkich wymuszeniach w parach obrotowych
    if os(wo(i,1),1) ~= 0
        Fq(rozF,os(wo(i,1),1)*3) = 1; % 2.42
    end
    if os(wo(i,1),2) ~= 0
        Fq(rozF,os(wo(i,1),2)*3) = -1; % 2.44
    end
    rozF = rozF + 1;
end

for i=1:size(wp, 1) % petla po wszystkich wymuszeniach w parach postepowych
    u = ps(wp(i,1),4:5)';
    sA = ps(wp(i,1),6:7)';
    [ri,~,Roti] = FromQ(q,ps(wp(i,1),1));
    [rj,~,Rotj] = FromQ(q,ps(wp(i,1),2));
    if ps(wp(i,1),1) ~= 0
        Fq(rozF,ps(wp(i,1),1)*3-2:ps(wp(i,1),1)*3-1) = -(Rotj * v)'; % 2.47
        Fq(rozF,ps(wp(i,1),1)*3) = -(Rotj * u)' * Om * Roti * sA; % 2.48
    end
    if ps(wp(i,1),2) ~= 0
        Fq(rozF,ps(wp(i,1),2)*3-2:ps(wp(i,1),2)*3-1) = (Rotj * v)'; % 2.49
        Fq(rozF,ps(wp(i,1),2)*3) = -(Rotj * u)' * Om * (rj - ri - Roti * sA); % 2.50
    end
    rozF = rozF + 1;
end
