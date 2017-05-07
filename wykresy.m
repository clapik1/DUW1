i=5;

subplot(3,1,1);
plot(t, q(3*i-2, :));
hold on;
grid minor;
plot(t, q(3*i-1, :));

subplot(3,1,2);
plot(t, dq(3*i-2, :));
hold on;
grid minor;
plot(t, dq(3*i-1, :));

subplot(3,1,3);
plot(t, ddq(3*i-2, :));
hold on;
grid minor;
plot(t, ddq(3*i-1, :));