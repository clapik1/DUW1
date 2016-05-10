function df=DWymuszenie(indeks,t)

switch indeks
    case 1
        df = -0.2/20 * cos(t/20);
    case 2
        df = -0.1/10 * cos(t/10);
    case 3
        df = -3/180*pi;
    case 4
        df = -2*t;
    case 5
        df= -1;
    case 6
        df= 0;    
    otherwise
        disp('fatal error in DWymuszenie.m')
end