%wczyt danych
c = fopen('ciala.txt', 'r');
p = fopen('pary.txt', 'r');
w = fopen('wymuszenia.txt', 'r');

%czytaj ciala
ilec=str2num(fgetl(c));
for i=1:ilec
    q0([3*(i-1)+1:3*i])=str2num(fgetl(c));
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
        qi=q0([3*i-2:3*i-1]);  fii=q0(3*i);
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
    %i j fi0 v punkt1 punkt2
    i=tmp(1); j=tmp(2);
    sA=[tmp(3); tmp(4)];
    sB=[tmp(5); tmp(6)];
    
    if(i==0)
        qi=[0;0]; fii=0;
    else
        qi=q0([3*i-2:3*i-1]);  fii=q0(3*i);
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
    
    v=[sA(2)-sB(2); sB(1)-sA(1)];
    v=v/norm(v);
    v=Rotj'*v;
    ps(k, 4:5)=v';
end

%czytaj wymuszenia obrotowe
ilewo=str2num(fgetl(w));
wo=zeros(ilewo, 2);
for k=1:ilewo
    wo(k,:)=str2num(fgetl(w));
end

%czytaj wymuszenia postepowe
ilewp=str2num(fgetl(w));
wp=zeros(ilewp, 2);
for k=1:ilewp
    wp(k,:)=str2num(fgetl(w));
end
