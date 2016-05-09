function f=Wymuszenie(indeks,t)

switch indeks
    case 1
        f = 10 + 5 * sin(10 * t);
    case 2
        f = t^2 + pi/2;
    otherwise
        f = 0;
        disp('fatal error in Wymuszenie.m')
end