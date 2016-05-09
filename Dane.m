% Dane.m � plik zawieraj�cy dane o wymiarach mechanizmu

a=1; b=5; c=4;  % Wymiary mechanizmu
sA0=[0;0]; sB1=[0;0];  % Para obrotowa 0-1
sA1=[a;0]; sB2=[0;0];  % Para obrotowa 1-2
sA2=[b;0]; sB3=[0;0];  % Para obrotowa 2-3
f0=0; v0=[0;1]; sA3=[0;0]; sB0=[0;c];  % Para post�powa 3-0

q0 = [0; a/2; 0; sqrt(b^2 - (c-a)^2)/2; (a+c)/2; 0; sqrt(b^2 - (c-a)^2); c; 0];
os = [0 1 0 0 0 -a/2; 1 2 0 a/2 (-sqrt(b^2-(c-a)^2)/2) (a-c)/2; 2 3 sqrt(b^2-(c-a)^2)/2 (c-a)/2 0 0];
ps = [3 0 0 0 1 0 0 0 c];
wo = [1 2];
wp = [];