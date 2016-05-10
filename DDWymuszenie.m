function ddf=DDWymuszenie(indeks,t)

switch indeks
    case 1
        ddf = - 0.2/20/20 * sin(t/20);
    case 2
        ddf = - 0.1/10/10 * sin(t/10);
    otherwise
        disp('fatal error in DDWymuszenie.m')
end