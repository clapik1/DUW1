function dq = Predkosc(q, os, ps, wo, wp, t)
%
%   Procedura do rozwi�zywania zadania o pr�dko�ci.
%   Zadanie o po�o�eniu musi by� rozwi�zane wcze�niej.
% Wej�cie:
%   q - wsp�rz�dne absolutne uk�adu wielocz�onowego,
%   t - chwila, dla kt�rej poszukiwane jest rozwi�zanie.
% Wyj�cie:
%   dq - obliczone pochodne wsp�rz�dnych absolutnych wzgl�dem czasu.
%

Ft = zeros(size(q, 1), 1);

ind = 2 * (size(os, 1) + size(ps, 1)) + 1;

for i=1:size(wo,1)
    Ft(ind, 1) = -DWymuszenie(wo(i, 2), t);
    ind = ind+1;
end

for i=1:size(wp,1)
    Ft(ind, 1) = -DWymuszenie(wp(i, 2), t);
    ind=ind+1;
end

Fq = Jakobian(q, os, ps, wo, wp);
dq = -Fq \ Ft;
