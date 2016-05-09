function F=Wiezy(q,t)
% F=Wiezy(q,t)
%   Procedura wsp�pracuj�ca z NewRaph.
%   S�u�y do obliczania warto�ci funkcji opisuj�cych wi�zy.
% Wej�cie:
%   q - wsp�rz�dne absolutne uk�adu wielocz�onowego,
%   t - aktualna chwila.
% Wyj�cie:
%   F - warto�ci funkcji.
%
%*************************************************************
%   Program stanowi za��cznik do ksi��ki:
%   Fr�czek J., Wojtyra M.: Kinematyka uk�ad�w wielocz�onowych.
%                           Metody obliczeniowe. WNT 2007.
%   Wersja 1.0
%*************************************************************

Dane; % Wczytanie danych o wymiarach mechanizmu i usytuowaniu par kinematycznych

% Przypisanie elementom wektora q czytelnych nazw
r1=q(1:2); fi1=q(3);   r2=q(4:5); fi2=q(6);   r3=q(7:8);   fi3=q(9);

% Obliczenie macierzy kosinus�w kierunkowych
Rot1=Rot(fi1);   Rot2=Rot(fi2);   Rot3=Rot(fi3);

% R�wnania wi�z�w
F(1:2,1)=sA0-(r1+Rot1*sB1);
F(3:4,1)=r1+Rot1*sA1-(r2+Rot2*sB2);
F(5:6,1)=r2+Rot2*sA2-(r3+Rot3*sB3);
F(7,1)=fi3-f0;
F(8,1)=v0'*(sB0-r3-Rot3*sA3);
F(9,1)=fi1-t^2-pi/2;
disp(F)

% wlasne - lepsze :D
F = ones(length(q0),1);

rozos = size(os);
rozF = 1;
for i=1:rozos(1) % petla po wszystkich parach obrotowych
    if os(i,1) == 0 % jezeli indeks mowi, ze chodzi o podstawe to przypisz zera do zmiennych
        r1 = [0; 0];
        fi1 = 0;
    else % w przeciwnym przypadku wczytaj zmienne z wektora q
        r1 = q(os(i,1)*3-2:os(i,1)*3-1);
        fi1 = q(os(i,1)*3);
    end
    Rot1 = Rot(fi1);
    if os(i,2) == 0 % to samo co wyzej dla drugiego czlonu pary
        r2 = [0; 0];
        fi2 = 0;
    else
        r2 = q(os(i,2)*3-2:os(i,2)*3-1);
        fi2 = q(os(i,2)*3);
    end
    Rot2 = Rot(fi2);
    F(rozF:rozF+1,1) = r1 + Rot1 * os(i,3:4)' - (r2 + Rot2 * os(i,5:6)'); % wzor 2.16 na pare obrotowa
    rozF = rozF + 2;
end

rozps = size(ps);
for i=1:rozps(1) % petla po wszystkich parach postepowych
    if ps(i,1) == 0
        r1 = [0; 0];
        fi1 = 0;
    else
        r1 = q(ps(i,1)*3-2:ps(i,1)*3-1);
        fi1 = q(ps(i,1)*3);
    end
    Rot1 = Rot(fi1);
    if ps(i,2) == 0
        r2 = [0; 0];
        fi2 = 0;
    else
        r2 = q(ps(i,2)*3-2:ps(i,2)*3-1);
        fi2 = q(ps(i,2)*3);
    end
    Rot2 = Rot(fi2);
    F(rozF,1) = fi1 - fi2 - ps(i,3); % wzor 2.17 na warunek kata w parze postepowej
    F(rozF+1,1) = (Rot2 * ps(i,4:5)')'*(r2 - r1 - Rot1 * ps(i,6:7)') + ps(i,4:5)*ps(i,8:9)'; % wzor 2.20 na brak ruchu w kierunku v
    rozF = rozF + 2;
end

rozwo = size(wo);
for i=1:rozwo % petla po wszystkich wymuszeniach w parach obrotowych
    if os(rozwo(i,1),1) == 0
        fi1 = 0;
    else
        fi1 = q(os(rozwo(i,1),1)*3);
    end
    if os(rozwo(i,1),2) == 0
        fi2 = 0;
    else
        fi2 = q(os(rozwo(i,1),2)*3);
    end
    F(rozF,1) = fi1 - fi2 - Wymuszenie(rozwo(i,2),t); % wzor 2.25 na wzajemna orientacje czlonow
    rozF = rozF + 1;
end

rozwp = size(wp);
for i=1:rozwp(1) % petla po wszystkich wymuszeniach w parach postepowych
    if os(rozwp(i,1),1) == 0
        r1 = [0; 0];
        fi1 = 0;
    else
        r1 = q(ps(rozwp(i,1),1)*3-2:ps(rozwp(i,1),1)*3-1);
        fi1 = q(ps(rozwp(i,1),1)*3);
    end
    Rot1 = Rot(fi1);
    if os(rozwp(i,1),2) == 0
        r2 = [0; 0];
        fi2 = 0;
    else
        r2 = q(ps(rozwp(i,1),2)*3-2:ps(rozwp(i,1),2)*3-1);
        fi2 = q(ps(rozwp(i,1),2)*3);
    end
    Rot2 = Rot(fi2);
    % wzor 2.26 na wymuszenie w parze postepowej
    F(rozF,1) = (Rot2 * ps(rozwp(i,1),4:5)')'*(r2 + Rot2 * ps(rozwp(i,1),8:9) - r1 - Rot1 * ps(rozwp(i,1),6:7)') - Wymuszenie(rozwp(i,2),t);
    rozF = rozF + 1;
end
