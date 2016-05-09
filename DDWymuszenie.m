function ddf=DDWymuszenie(indeks,t)

switch indeks
    case 1
        ddf = -500 * sin(10 * t);
    otherwise
        ddf = 0;
        disp('fatal error in DDWymuszenie.m')
end