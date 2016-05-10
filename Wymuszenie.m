function f=Wymuszenie(indeks,t)

switch indeks
    case 1
        f = sqrt(0.2^2 + 0.6^2) + 0.2 * sin(t/20);
    case 2
        f = sqrt(0.2^2 + 0.8^2) + 0.1 * sin(t/10);
    otherwise
        f = 0;
        disp('fatal error in Wymuszenie.m')
end