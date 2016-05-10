function [r, fi, Roti] = FromQ(q, i)

if i == 0 % jezeli indeks mowi, ze chodzi o podstawe to przypisz zera do zmiennych
    r = [0; 0];
    fi = 0;
else % w przeciwnym przypadku wczytaj zmienne z wektora q
    r = q(i*3-2:i*3-1);
    fi = q(i*3);
end
Roti = Rot(fi);