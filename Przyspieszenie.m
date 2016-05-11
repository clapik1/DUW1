function d2q = Przyspieszenie(dq, q, os, ps, wo, wp, t)
% d2q=Przyspieszenie(dq,q,t)
%   Procedura do rozwi¹zywania zadania o przyspieszeniu.
%   Zadania o po³o¿eniu i prêdkoci musz¹ byæ rozwi¹zane wczeniej.
% Wejcie:
%   dq - pochodne wspó³rzêdnych absolutnych wzglêdem czasu,
%   q  - wspó³rzêdne absolutne uk³adu wielocz³onowego,
%   t  - chwila, dla której poszukiwane jest rozwi¹zanie.
% Wyjcie:
%   d2q - obliczone drugie pochodne wspó³rzêdnych absolutnych wzglêdem czasu.
%

Om = [0 -1; 1 0]; % Stala macierz

gam = zeros(length(q), 1);

ind=1;

for k=1:size(os,1) %petla po wszystkich parach obrotowych
    i = os(k, 1);
    j = os(k, 2);
    
    [~,~,Roti] = FromQ(q, i);
    [~,~,Rotj] = FromQ(q, j);
    
    sA = os(k, 3:4)';
    sB = os(k, 5:6)';
    
    [~,dfii,~]=FromQ(dq,i);
    [~,dfij,~]=FromQ(dq,j);
    
    gam(ind:ind+1, 1) = Roti * sA * dfii^2 - Rotj * sB * dfij^2; % wzor (2.40)
    
    ind = ind+2;
end

for k=1:size(ps,1) %petla po wszystkich parach postepowych
    gam(ind, 1) = 0; % wzor (2.46)
    ind=ind+1;
    
    i = ps(k, 1);
    j = ps(k, 2);
    
    u = ps(k, 4:5)';
    v = [-u(2); u(1)];
    sA = ps(k, 6:7)';
    
    [ri,~,Roti] = FromQ(q, i);
    [rj,~,Rotj] = FromQ(q, j);
    
    [dri,dfii,~]=FromQ(dq,i);
    [drj,dfij,~]=FromQ(dq,j);
    
    gam(ind, 1) = (Rotj * v)'*(2 * Om * (drj - dri) * dfij + (rj - ri) * dfij^2 - Roti * sA * (dfij - dfii)^2 ); % ... wzor (2.57)
    
    ind = ind + 1;
end

for k=1:size(wo,1) %petla po wszystkich wymuszeniach obrotowych
    gam(ind,1) = -DDWymuszenie(wo(k, 2), t);
    ind=ind+1;
end

for k=1:size(wp,1) %petla po wszystkich wymuszeniach postepowych    
    i = ps(wp(k, 1), 1);
    j = ps(wp(k, 1), 2);
    
    u = ps(wp(k, 1), 4:5)';
    sA = ps(wp(k, 1), 6:7)';
    
    [ri,~,Roti] = FromQ(q, i);
    [rj,~,Rotj] = FromQ(q, j);
    
    [dri,dfii,~]=FromQ(dq,i);
    [drj,dfij,~]=FromQ(dq,j);
    
    gam(ind, 1) = (Rotj * u)'*(2 * Om * (drj - dri) * dfij + (rj - ri) * dfij^2 - Roti * sA * (dfij - dfii)^2 ) + DDWymuszenie(wp(k, 2), t); % ... wzor (2.57) + pochodna wymuszenia
    
    ind = ind + 1;
end

% Obliczenie macierzy uk³adu równañ
Fq = Jakobian(q, os, ps, wo, wp);

% Obliczenie przyspieszenia (rozwi¹zanie uk³adu równañ)
d2q = Fq \ gam;
