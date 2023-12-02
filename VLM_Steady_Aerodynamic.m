clear; clc;

fprintf('Calculation Start\n');

%% parameters
m_x = 20;                      % number of chord-wise panels on the wing
n_x = 5*m_x + m_x;             % total number of chord-wise panels
n_y = 20;                      % total number of span-wise panels 

AR = 1.27;
chord = 1;                   % chord of the wing
span = chord * AR;           % span of the wing
  
u     = 10;                     % velocity
RF    = 0.992;                 % relaxation factor
ro    = 1.225;              % density
n_t   = 1000;                  % number of time steps
alpha = 4 *pi/180;                     % AOA

%% -------- %%
m = m_x * (n_y);
n = n_x * (n_y);

sum_gammas = 1:n;
sum_gammas = reshape(sum_gammas,[n_y,n_x]).';
 

A_plate = chord * span;
dx = chord/m_x;
dy = span/n_y;

dt = dx/u;
time = 0:dt:(n_t-1)*dt;
s_time = time*u/chord;


%% index
N = [];
for i = 1:n_x
    for j = 1:n_y
        N = [N;i,j];
    end
end

%% get Xs & Ys
for i = 1 : n_x
    for j = 1 : n_y
        x_i(i,j) = ((i - 1)*dx + 0.75*dx);
        x_j(i,j) = ((i - 1)*dx + 0.25*dx);
        y_i(i,j) = j * dy - 0.5 * dy - span/2;
        y_jb(i,j) = j * dy - span/2;
        y_ja(i,j) = (j-1) * dy - span/2;
    end
end

x_interval = x_i(1:m_x);
y_interval = y_i(1,:);

%% get K
for i = 1:m
    for j = 1:n
        K_a(i,j) = 1 + sqrt((x_i(N(i,1),N(i,2)) - x_j(N(j,1),N(j,2)))^2 + (y_i(N(i,1),N(i,2)) - y_ja(N(j,1),N(j,2)))^2) / (x_i(N(i,1),N(i,2)) - x_j(N(j,1),N(j,2)));
        K_b(i,j) = 1 + sqrt((x_i(N(i,1),N(i,2)) - x_j(N(j,1),N(j,2)))^2 + (y_i(N(i,1),N(i,2)) - y_jb(N(j,1),N(j,2)))^2) / (x_i(N(i,1),N(i,2)) - x_j(N(j,1),N(j,2)));
        K(i,j) = -1/(4*pi*(y_i(N(i,1),N(i,2)) - y_ja(N(j,1),N(j,2)))) * K_a(i,j) + 1/(4*pi*(y_i(N(i,1),N(i,2)) - y_jb(N(j,1),N(j,2)))) * K_b(i,j);
    end
end

%% Local Slope
ind = 1;
% ---- Airfoil Camber Profile ---- for symmetric airfoils: z_c = zeros(1,m_x)
while x_interval(ind) <0.4
    z_c(ind) = 0.125*(0.8*x_interval(ind) - (x_interval(ind))^2);
    ind = ind+1;
end
for i = ind:m_x
    z_c(i) = 0.0556* (0.2+0.8*x_interval(i) - (x_interval(i))^2);
end
% ---- ---- ---- ---- ---- ---- ----

a_local = atan(diff([0,z_c]) ./ diff([0,x_interval]));

%% get w
a_net = alpha + repelem(a_local, n_y);

for t = 1:n_t
    w_m = -1*a_net*u;
    w(:,t) = [w_m'; zeros(n-m, 1)];
end

%% get A
I = eye(n_y);
A_1 = [];
for i = 1:m_x
    A_1 = [A_1 I];
end
A = [K;
    A_1, I, zeros(n_y, n-m-2*n_y), zeros(n_y, n_y);
    zeros(n-m-2*n_y, m), zeros(n-m-2*n_y, n_y), eye(n-m-2*n_y), zeros(n-m-2*n_y, n_y);
    zeros(n_y, m), zeros(n_y, n_y), zeros(n_y, n-m-2*n_y), eye(n_y)];

%% get B
B_1 = -1*A_1;
B_2 = [zeros(n-m-2*n_y, n_y);
      -1*RF*eye(n_y)];

B = [zeros(m, n);
    B_1, zeros(n_y), zeros(n_y, n-m-2*n_y), zeros(n_y, n_y);
    zeros(n-m-n_y, m), -1*eye(n-m-n_y), B_2];

[L,U] = lu(A);


%% get Gamma
gamma_old = zeros(n, n_t);
gamma = zeros(n, n_t);
for t = 0 : n_t
    if t == 0
        for i = 1 : m
            w_m0(i) = -1 * u;
        end
        w0 = [w_m0'; zeros(n-m, 1)];
        RHS = w0;
    else
        RHS = w(:,t) - B * gamma(:,t);
    end
    Y = L\RHS;
    gamma(:,t+1) = U\Y;
end



%% Get P & L
for t = 0 : n_t
    count = 0;
    if t == 0
        for ii = 1:m
            if mod(ii,n_y) == 0
                count = count + 1;
            end
            if mod(ii,n_y) == 0
                P_nd(ii,t+1) = 1/u/dx*((gamma(ii,t+1)+0)/2 + (sum(gamma(sum_gammas(1:count,n_y),t+1) - 0)));
                L_nd(ii,t+1) = ro*u*gamma(ii,t+1)*dy + ro*dx*dy*(sum(gamma(sum_gammas(1:count,n_y),t+1) - 0));
            else
                P_nd(ii,t+1) = 1/u/dx*((gamma(ii,t+1)+0)/2 + (sum(gamma(sum_gammas(1:count,mod(ii,n_y)),t+1) - 0)));
                L_nd(ii,t+1) = ro*u*gamma(ii,t+1)*dy + ro*dx*dy*(sum(gamma(sum_gammas(1:count,mod(ii,n_y)),t+1) - 0));
            end
        end
    else
        for ii = 1:m

            if mod(ii,n_y) == 0
                count = count + 1;
            end

            if mod(ii,n_y) == 0
                P_nd(ii,t+1) = 1/u/dx*((gamma(ii,t+1)+gamma(ii,t))/2 + (sum(gamma(sum_gammas(1:count,n_y),t+1) - gamma(sum_gammas(1:count,n_y),t))));
                L_nd(ii,t+1) = ro*u*gamma(ii,t+1)*dy + ro*dx*dy*(sum(gamma(sum_gammas(1:count,n_y),t+1) - gamma(sum_gammas(1:count,n_y),t)) + 3/4*(gamma(ii,t+1) - gamma(ii,t)));
            else
                P_nd(ii,t+1) = 1/u/dx*((gamma(ii,t+1)+gamma(ii,t))/2 + (sum(gamma(sum_gammas(1:count,mod(ii,n_y)),t+1) - gamma(sum_gammas(1:count,mod(ii,n_y)),t))));
                L_nd(ii,t+1) = ro*u*gamma(ii,t+1)*dy + ro*dx*dy*(sum(gamma(sum_gammas(1:count,mod(ii,n_y)),t+1) - gamma(sum_gammas(1:count,mod(ii,n_y)),t)) + 3/4*(gamma(ii,t+1) - gamma(ii,t)));
            end
        end
    end
end

L_nd(:,1) = [];
L_ND = sum(L_nd(1:m),1);
P_ND = -sum(P_nd,1);
C_L = sum(L_nd,1)/(0.5*ro*u^2*A_plate);


P_1 = P_nd(:,end);

P_grid = [];
for i = 1 : m_x
    P_grid = [P_grid; P_1((i-1)*n_y+1 : i*n_y)'];
end
P_grid = P_grid';
P_chord = trapz(y_interval, P_grid);
P_span = trapz(x_interval, P_grid');

%% Center of Pressure
x_interval = x_i(1:m_x);
y_interval = y_i(1,:);
x_cp = trapz(x_interval, x_interval .* P_chord) / trapz(x_interval, P_chord);
y_cp = trapz(y_interval, y_interval .* P_span) / trapz(y_interval, P_span);

%% Plot
figure();
surf(x_interval, y_interval, P_grid*2); grid on;