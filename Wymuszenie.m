function f=Wymuszenie(indeks,t)

switch indeks
    case 1
        f = sqrt(0.2^2 + 0.6^2) + 0.2 * sin(t/20);
        %f = 0.2 * sin(t/20);
    case 2
        f = sqrt(0.2^2 + 0.8^2) + 0.1 * sin(t/10);
        %f = 0.1 * sin(t/10);
    case 3
        f = -3/180*pi*t;
    case 4
        f = -t^2 - pi/2;
    otherwise
        disp('fatal error in Wymuszenie.m')
end