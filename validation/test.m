

n_t = 200;
u = ones(n_t, 1);
t = linspace(0, 10, n_t);

dt = t(2) - t(1);
fs = 1/dt;

tf1 = tf([1 1], [1 1 1]);
tf2 = tf([1], [1 1 1]);

y1 = lsim(tf1, u, t);
y2 = lsim(tf2, u, t);

figure()
hold on
plot(t, y1)
plot(t, y2)
plot(t, u, "--k")
legend("1", "2", "input")
hold off
