function df=DWymuszenie(indeks,t)

switch indeks
    case 1
        df = 50 * cos(10 * t);
        
    otherwise
        df = 0;
        disp('fatal error in DWymuszenie.m')
end