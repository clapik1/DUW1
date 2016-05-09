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
    if os(i,1) == 0 % jezeli indeks mowi, ze chodzi o podstawe to przypisz zera do zmiennych
        ri = [0; 0];
        fii = 0;
    else % w przeciwnym przypadku wczytaj zmienne z wektora q
        ri = q(os(i,1)*3-2:os(i,1)*3-1);
        fii = q(os(i,1)*3);
    end
    Roti = Rot(fii);
    if os(i,2) == 0 % to samo co wyzej dla drugiego czlonu pary
        rj = [0; 0];
        fij = 0;
    else
        rj = q(os(i,2)*3-2:os(i,2)*3-1);
        fij = q(os(i,2)*3);
    end
    Rotj = Rot(fij);
    sA = os(i,3:4)';
    sB = os(i,5:6)';
    F(rozF:rozF+1,1) = ri + Roti * sA - (rj + Rotj * sB); % wzor 2.16 na pare obrotowa
    rozF = rozF + 2;
end

rozps = size(ps);
for i=1:rozps(1) % petla po wszystkich parach postepowych
    if ps(i,1) == 0
        ri = [0; 0];
        fii = 0;
    else
        ri = q(ps(i,1)*3-2:ps(i,1)*3-1);
        fii = q(ps(i,1)*3);
    end
    Roti = Rot(fii);
    if ps(i,2) == 0
        rj = [0; 0];
        fij = 0;
    else
        rj = q(ps(i,2)*3-2:ps(i,2)*3-1);
        fij = q(ps(i,2)*3);
    end
    Rotj = Rot(fij);
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
    if os(rozwo(i,1),1) == 0
        fii = 0;
    else
        fii = q(os(rozwo(i,1),1)*3);
    end
    if os(rozwo(i,1),2) == 0
        fij = 0;
    else
        fij = q(os(rozwo(i,1),2)*3);
    end
    F(rozF,1) = fii - fij - Wymuszenie(rozwo(i,2),t); % wzor 2.25 na wzajemna orientacje czlonow
    rozF = rozF + 1;
end

rozwp = size(wp);
for i=1:rozwp(1) % petla po wszystkich wymuszeniach w parach postepowych
    if os(rozwp(i,1),1) == 0
        ri = [0; 0];
        fii = 0;
    else
        ri = q(ps(rozwp(i,1),1)*3-2:ps(rozwp(i,1),1)*3-1);
        fii = q(ps(rozwp(i,1),1)*3);
    end
    Roti = Rot(fii);
    if os(rozwp(i,1),2) == 0
        rj = [0; 0];
        fij = 0;
    else
        rj = q(ps(rozwp(i,1),2)*3-2:ps(rozwp(i,1),2)*3-1);
        fij = q(ps(rozwp(i,1),2)*3);
    end
    Rotj = Rot(fij);
    u = ps(rozwp(i,1),4:5)';
    sA = ps(rozwp(i,1),6:7)';
    sB = ps(rozwp(i,1),8:9);
    % wzor 2.26 na wymuszenie w parze postepowej
    F(rozF,1) = (Rotj * u)'*(rj + Rotj * sB - ri - Roti * sA) - Wymuszenie(rozwp(i,2),t);
    rozF = rozF + 1;
end
