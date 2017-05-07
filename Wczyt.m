function [q0, m, os, ps, ws, wst, de] = Wczyt()

%wczyt danych
c = fopen('dane/ciala.txt', 'r');
p = fopen('dane/pary.txt', 'r');
w = fopen('dane/wymuszenia.txt', 'r');
d = fopen('dane/doexportu.txt', 'r');

%czytaj ciala
ilec=str2num(fgetl(c));
for i=1:ilec
    tmp = str2num(fgetl(c));
    q0(3*(i-1)+1:3*i)=tmp(1:3);
    m(i, 1:2)=tmp(4:5);
end
q0=q0';

%czytaj pary obrotowe
ileo=str2num(fgetl(p));
os=zeros(ileo, 6);
for k=1:ileo
    tmp=str2num(fgetl(p));
    i=tmp(1); j=tmp(2);
    s=[tmp(3); tmp(4)];
    
    if(i==0)
        qi=[0;0]; fii=0;
    else
        qi=q0(3*i-2:3*i-1);  fii=q0(3*i);
    end
       
    if(j==0)
        qj=[0;0]; fij=0;
    else
        qj=q0([3*j-2:3*j-1]); fij=q0(3*j);
    end
    
    os(k, 1:2)=[i j];
    
    Roti=Rot(fii);
    os(k, 3:4)=Roti'*(s-qi);
    
    Rotj=Rot(fij);
    os(k, 5:6)=Rotj'*(s-qj);
end

%czytaj pary postepowe
ilep=str2num(fgetl(p));
ps=zeros(ilep, 9);
for k=1:ilep
    tmp=str2num(fgetl(p));
    %i j fi0 _____u______ punkt1 punkt2
    i=tmp(1); j=tmp(2);
    sA=[tmp(3); tmp(4)];
    sB=[tmp(5); tmp(6)];
    
    if(i==0)
        qi=[0;0]; fii=0;
    else
        qi=q0(3*i-2:3*i-1);  fii=q0(3*i);
    end
       
    if(j==0)
        qj=[0;0]; fij=0;
    else
        qj=q0([3*j-2:3*j-1]); fij=q0(3*j);
    end
    
    ps(k, 1:2)=[i j];
    ps(k, 3)=fij-fii;
    
    Roti=Rot(fii);
    ps(k, 6:7)=Roti'*(sA-qi);
    
    Rotj=Rot(fij);
    ps(k, 8:9)=Rotj'*(sB-qj);
    
    u=[sB(1)-sA(1); sB(2)-sA(2)];
    u=u/norm(u);
    u=Rotj'*u;
    ps(k, 4:5)=u';
end

%czytaj przylozone sily
ilews=str2num(fgetl(w));
ws=zeros(ilews, 5); %i x_i y_i F_x F_y
for k=1:ilews
    tmp=str2num(fgetl(w)); %i x_w_glob y_w_glob F fi_F
    
    i=tmp(1);
    s=[tmp(2); tmp(3)];
    F=[tmp(4) 0]'; fi_F=deg2rad(tmp(5));
    
    ws(k, 1)=i;
    
    qi=q0(3*i-2:3*i-1);  fii=q0(3*i);
    
    Roti=Rot(fii);
    ws(k, 2:3)=Roti'*(s-qi);
    
    ws(k, 4:5)=Rot(fi_F)*F;
end

%czytaj elementy sprezysto tlumiace
ilewst=str2num(fgetl(w));
wst=zeros(ilewst, 9); %i j x_i y_i x_j y_j spr l0 tlum
for k=1:ilewst
    tmp=str2num(fgetl(w));
    
    i=tmp(1); j=tmp(2);
    sA=[tmp(3); tmp(4)];
    sB=[tmp(5); tmp(6)];
    spr=tmp(7); tlum=tmp(8);
    
    if(i==0)
        qi=[0;0]; fii=0;
    else
        qi=q0(3*i-2:3*i-1);  fii=q0(3*i);
    end
       
    if(j==0)
        qj=[0;0]; fij=0;
    else
        qj=q0([3*j-2:3*j-1]); fij=q0(3*j);
    end
    
    wst(k, 1:2)=[i j];
    
    Roti=Rot(fii);
    wst(k, 3:4)=Roti'*(sA-qi);
    
    Rotj=Rot(fij);
    wst(k, 5:6)=Rotj'*(sB-qj);
    
    wst(k, 7:9)=[spr norm(sB-sA) tlum]';
end

%czytaj co eksportowac
ilepun=str2num(fgetl(d));
de=zeros(ilepun, 3);
for k=1:ilepun
    tmp=str2num(fgetl(d));
   
    i=tmp(1);
    de(k, 1)=i;
    
    s=[tmp(2); tmp(3)];
    
    qi=q0(3*i-2:3*i-1);
    fii=q0(3*i);
    
    Roti=Rot(fii);
    de(k, 2:3)=Roti'*(s-qi);
end
