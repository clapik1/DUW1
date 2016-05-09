function Fq=Jakobian(q)
% Fq=Jakobian(q)
%   Procedura wspolpracujaca z NewRaph.
%   Sluzy do obliczania macierzy Jacobiego rownan wiezow.
% Wejscie:
%   q - wspolrzedne absolutne ukladu wieloczlonowego.
% Wyjscie:
%   Fq - obliczona macierz Jacobiego.
%

Dane; % Wczytanie danych o wymiarach mechanizmu i usytuowaniu par kinematycznych

Om=[0 -1;1 0];  % Sta�a macierz

% Przypisanie elementom wektora q czytelnych nazw
r1=q(1:2); fi1=q(3);   r2=q(4:5); fi2=q(6);  r3=q(7:8); fi3=q(9);

% Obliczenie macierzy kosinus�w kierunkowych
Rot1=Rot(fi1);  Rot2=Rot(fi2);  Rot3=Rot(fi3);

% Macierz Jaocobiego - pocz�tkowo zerowa
Fq=zeros(9,9);
% Macierz Jaocobiego - wype�nianie niezerowych element�w
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

% wlasne!

rozos = size(os);
rozF = 1;
for i=1:rozos(1) % petla po wszystkich parach obrotowych
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

rozps = size(ps);
for i=1:rozps(1) % petla po wszystkich parach postepowych
    v = ps(i,4:5)';
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

rozwo = size(wo);
for i=1:rozwo % petla po wszystkich wymuszeniach w parach obrotowych
    if os(rozwo(i,1),1) ~= 0
        Fq(rozF,os(rozwo(i,1),1)*3) = 1; % 2.42
    end
    if os(rozwo(i,1),2) ~= 0
        Fq(rozF,os(rozwo(i,1),2)*3) = -1; % 2.44
    end
    rozF = rozF + 1;
end

rozwp = size(wp);
for i=1:rozwp(1) % petla po wszystkich wymuszeniach w parach postepowych
    v = ps(rozwp(i,1),4:5)';
    sA = ps(rozwp(i,1),6:7)';
    [ri,~,Roti] = FromQ(q,ps(rozwp(i,1),1));
    [rj,~,Rotj] = FromQ(q,ps(rozwp(i,1),2));
    if ps(rozwp(i,1),1) ~= 0
        Fq(rozF,ps(rozwp(i,1),1)*3-2:ps(rozwp(i,1),1)*3-1) = -(Rotj * v)'; % 2.47
        Fq(rozF,ps(rozwp(i,1),1)*3) = -(Rotj * v)' * Om * Roti * sA; % 2.48
    end
    if ps(rozwp(i,1),2) ~= 0
        Fq(rozF,ps(rozwp(i,1),2)*3-2:ps(rozwp(i,1),2)*3-1) = (Rotj * v)'; % 2.49
        Fq(rozF,ps(rozwp(i,1),2)*3) = -(Rotj * v)' * Om * (rj - ri - Roti * sA); % 2.50
    end
    rozF = rozF + 1;
end
