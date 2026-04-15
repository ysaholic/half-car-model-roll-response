%% --- Parameters ---
ms = 750 ;     %% Sprung mass                               [kg]
mu = 35 ;      %% Unsprung mass                             [kg]
I = 320 ;      %% Mass moment of inertia                    [kgm2]
kt = 200000 ;     %% spring const. wheel & tire             [N/m]
ks = 20000 ;      %% spring const. suspension               [N/m]
cs = 1500 ;    %% damping const. suspension                 [Ns/m]
l = 0.3 ;      %% distance between COG & suspension         [Ns/m]

mur = mu;  mul = mu;
ksr = ks;  ksl = ks;
csr = cs;  csl = cs;
ktr = kt;  ktl = kt;
lr  = l ;  ll  = l ;

%% --- State-Space Model ---

% --- A matrix: [A1 A2; A3 A4] ---
A1 = zeros(4,4);
A2 = eye(4);

% A3: displacement (stiffness) coupling terms
A3 = [
    -(ksr+ktr)/mur,        0,              ksr/mur,          ksr*lr/mur;
         0,          -(ksl+ktl)/mul,        ksl/mul,         -ksl*ll/mul;
    ksr/ms,               ksl/ms,       -(ksr+ksl)/ms,   -(ksr*lr-ksl*ll)/ms;
    ksr*lr/I,            -ksl*ll/I,   -(ksr*lr-ksl*ll)/I, -(ksr*lr^2+ksl*ll^2)/I
];

% A4: velocity (damping) coupling terms
A4 = [
    -csr/mur,       0,          csr/mur,         csr*lr/mur;
        0,       -csl/mul,       csl/mul,        -csl*ll/mul;
    csr/ms,       csl/ms,    -(csr+csl)/ms,  -(csr*lr-csl*ll)/ms;
    csr*lr/I,   -csl*ll/I, -(csr*lr-csl*ll)/I, -(csr*lr^2+csl*ll^2)/I
];

A = [A1, A2;
     A3, A4];

% --- B matrix: control input u = [ur; ul] = [0; 0] (passive) ---
B = [
     0,       0;
     0,       0;
     0,       0;
     0,       0;
    -1/mur,   0;
     0,      -1/mul;
     1/ms,    1/ms;
     lr/I,   -ll/I
];

% --- G matrix: road disturbance input z_road = [zrr; zrl] ---
G = [
     0,       0;
     0,       0;
     0,       0;
     0,       0;
    ktr/mur,  0;
     0,      ktl/mul;
     0,       0;
     0,       0
];

% --- C matrix: output, 8 state  ---
C = eye(8);
D = zeros(8,2);
H = zeros(8,2);

% Display matriks
fprintf('=== State-Space Matrices ===\n');
fprintf('\nA matrix (8x8):\n'); disp(A);
fprintf('B matrix (8x2) [control input]:\n'); disp(B);
fprintf('G matrix (8x2) [road disturbance input]:\n'); disp(G);

% Cek eigenvalue (cek kestabilan)
ev = eig(A);
fprintf('Eigenvalues of A:\n');
disp(ev);
if all(real(ev) < 0)
    fprintf('-> STABLE (all eigenvalues have negative real parts)\n\n');
else
    fprintf('-> WARNING: System may be UNSTABLE\n\n');
end
%% --- RK4 Simulation ---
t_end = 10;         % Total simulation time [s]
dt    = 0.001;      % Time step             [s]
t     = 0:dt:t_end;
N     = length(t);

% --- Build road profile vectors ---
zrr = zeros(1,N);
zrl = zeros(1,N);
for i = 1:N
    ti = t(i);
    if ti >= t1_r && ti <= (t1_r + T_bump)
        zrr(i) = (a_r/2) * (1 - cos(omega * (ti - t1_r)));
    end
    if ti >= t1_l && ti <= (t1_l + T_bump)
        zrl(i) = (a_l/2) * (1 - cos(omega * (ti - t1_l)));
    end
end

% --- RK4 time integration ---
x      = zeros(8, N);
x(:,1) = zeros(8,1);     % Initial condition: system at rest
u_ctrl = zeros(2,1);     % Passive: no control force

for i = 1:N-1
    z_i = [zrr(i); zrl(i)];

    f  = @(xi) A*xi + B*u_ctrl + G*z_i;

    k1 = f(x(:,i));
    k2 = f(x(:,i) + (dt/2)*k1);
    k3 = f(x(:,i) + (dt/2)*k2);
    k4 = f(x(:,i) +  dt   *k3);

    x(:,i+1) = x(:,i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

%%  ---  ROAD INPUT PROFILE ---
%  Half-sine bump
%    z_r(t) = (a/2) * (1 - cos(2*pi*V/b * (t-t1)))  for t1 <= t <= t1 + b/V
%             0

a_r  = 0.05;   % Right wheel bump amplitude  [m]
a_l  = 0.08;   % Left wheel bump amplitude   [m] 
b = 0.5        % Width of half-sine bump     [m]
V    = 10;     % Vehicle speed               [m/s]
t1_r = 1.0;    % Bump hit time, right wheel  [s]
t1_l = 1.3;    % Bump hit time, left wheel   [s]

T_bump = b / V;    % Duration of bump excitation [s]
omega  = 2*pi*V/b; % Bump angular frequency      [rad/s]
%% -------------------------------------------------------------------------

t_end = 10;         % Total simulation time [s]
dt    = 0.001;      % Time step             [s]
t     = 0:dt:t_end;
N     = length(t);

% --- Vektor profil jalan ---
zrr = zeros(1,N);
zrl = zeros(1,N);
for i = 1:N
    ti = t(i);
    if ti >= t1_r && ti <= (t1_r + T_bump)
        zrr(i) = (a_r/2) * (1 - cos(omega * (ti - t1_r)));
    end
    if ti >= t1_l && ti <= (t1_l + T_bump)
        zrl(i) = (a_l/2) * (1 - cos(omega * (ti - t1_l)));
    end
end

% --- RK4 time integration ---
x      = zeros(8, N);
x(:,1) = zeros(8,1);     % Initial condition: system at rest
u_ctrl = zeros(2,1);     % Passive: no control force

for i = 1:N-1
    z_i = [zrr(i); zrl(i)];

    f  = @(xi) A*xi + B*u_ctrl + G*z_i;

    k1 = f(x(:,i));
    k2 = f(x(:,i) + (dt/2)*k1);
    k3 = f(x(:,i) + (dt/2)*k2);
    k4 = f(x(:,i) +  dt   *k3);

    x(:,i+1) = x(:,i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end
%% --- OUTPUT ---
 
zur       = x(1,:);             % Right unsprung mass vertical disp.  [m]
zul       = x(2,:);             % Left  unsprung mass vertical disp.  [m]
zs        = x(3,:);             % Sprung mass vertical disp. (bounce) [m]
theta     = x(4,:);             % Roll angle                          [rad]
theta_deg = rad2deg(theta);     % Roll angle                          [deg]
zdot_s    = x(7,:);             % Sprung mass vertical velocity       [m/s]

% Suspension attachment points on sprung mass
zsr = zs + lr*theta;   % Right suspension contact point [m]
zsl = zs - ll*theta;   % Left  suspension contact point [m]

% Suspension deflection (working space indicator)
susp_def_r = zur - zsr;
susp_def_l = zul - zsl;

% Tire deflection (road holding indicator)
tire_def_r = zur - zrr;
tire_def_l = zul - zrl;
%% --- PERFORMANCE METRICS ---
fprintf('=== Performance Metrics ===\n');
fprintf('Max sprung mass displacement  : %8.4f mm\n',  max(abs(zs))      *1000);
fprintf('Max roll angle                : %8.4f deg\n', max(abs(theta_deg))     );
fprintf('Max suspension travel (right) : %8.4f mm\n',  max(abs(susp_def_r))*1000);
fprintf('Max suspension travel (left)  : %8.4f mm\n',  max(abs(susp_def_l))*1000);
fprintf('Max tire deflection  (right)  : %8.4f mm\n',  max(abs(tire_def_r))*1000);
fprintf('Max tire deflection  (left)   : %8.4f mm\n',  max(abs(tire_def_l))*1000);
%% --- PLOT ---
% --- Figure 1: Road Input ---
figure('Name','Road Input Profile','NumberTitle','off','Position',[100 100 700 300]);
plot(t, zrr*1000, 'b',  'LineWidth', 1.8); hold on;
plot(t, zrl*1000, 'r--','LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Road Profile z_r (mm)');
title('Road Disturbance Input (Different Roughness Left vs. Right)');
legend('Right wheel (z_{rr})','Left wheel (z_{rl})','Location','best');
grid on;

% --- Figure 2: Roll Response ---
figure('Name','Roll Response','NumberTitle','off','Position',[100 450 900 600]);

subplot(3,2,1)
plot(t, zs*1000, 'b','LineWidth',1.5);
xlabel('Time (s)'); ylabel('z_s (mm)');
title('Sprung Mass Vertical Displacement (Bounce)');
grid on;

subplot(3,2,2)
plot(t, theta_deg, 'r','LineWidth',1.5);
xlabel('Time (s)'); ylabel('\theta (deg)');
title('Roll Angle \theta (Sprung Mass)');
grid on;

subplot(3,2,3)
plot(t, zur*1000, 'b', 'LineWidth',1.5); hold on;
plot(t, zul*1000, 'r--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('z_u (mm)');
title('Unsprung Mass Displacement');
legend('Right (z_{ur})','Left (z_{ul})','Location','best');
grid on;

subplot(3,2,4)
plot(t, susp_def_r*1000, 'b', 'LineWidth',1.5); hold on;
plot(t, susp_def_l*1000, 'r--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Deflection (mm)');
title('Suspension Working Space (z_u - z_{sr/sl})');
legend('Right','Left','Location','best');
grid on;

subplot(3,2,5)
plot(t, tire_def_r*1000, 'b', 'LineWidth',1.5); hold on;
plot(t, tire_def_l*1000, 'r--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Tire Defl. (mm)');
title('Tire Deflection (Road Holding)');
legend('Right (z_{ur}-z_{rr})','Left (z_{ul}-z_{rl})','Location','best');
grid on;

subplot(3,2,6)
plot(t, zdot_s*1000, 'k','LineWidth',1.5);
xlabel('Time (s)'); ylabel('dz_s/dt (mm/s)');
title('Sprung Mass Vertical Velocity');
grid on;

sgtitle('Half-Car Passive Suspension — Roll Response (u_l = u_r = 0)');




