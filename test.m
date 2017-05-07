[q0, m, os, ps, ws, wst, de] = Wczyt();
dq0=zeros(size(q0,1),1);

%[F,Fq,G,M,Q]=Macierze(q0,dq0,os,ps,ws,wst,m);

tK=5.;
[t, q, dq, ddq] = Calkuj_Baumgarte(q0,tK,os,ps,ws,wst,m);
lroz = 0;

EXPQ = zeros(2, length(t));
EXPDQ = zeros(2, length(t));
EXPDDQ = zeros(2, length(t));
Om = [0 -1; 1 0];

for tt = 1:length(t)
    lroz = lroz + 1;
    qt = q(:,lroz);
    dqt = dq(:,lroz);
    ddqt = ddq(:,lroz);
    
    for i = 1:size(de, 1)
        sA = de(i,2:3)';
        [ri,~,Roti] = FromQ(qt, de(i, 1));
        [dri,dfii,~] = FromQ(dqt, de(i, 1));
        [ddri,ddfii,~] = FromQ(ddqt, de(i, 1));
        EXPQ(:, lroz) = ri + Roti * sA;
        EXPDQ(:, lroz) = dri + Om * Roti * sA * dfii;
        EXPDDQ(:, lroz) = ddri + Om * Roti * sA * ddfii - Roti * sA * dfii^2;
    end
end