function df=DWymuszenie(indeks,t)

switch indeks
    case 1
        df = 0.2/20 * cos(t/20);
    case 2
        df = 0.1/10 * cos(t/10);
    otherwise
        disp('fatal error in DWymuszenie.m')
end