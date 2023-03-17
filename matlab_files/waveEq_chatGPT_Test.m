clear all
close all
clc

% Define spatial and temporal discretization
L = 1; % length of the domain
T = 1; % total time
Nx = 100; % number of spatial points
Nt = 1000; % number of time steps

% Define the spatial and time step sizes
dx = L/Nx;
dt = T/Nt;

% Define the wave speed and the coefficient used in the finite difference method
c = 1;
coeff = c^2*dt^2/dx^2;

% Initialize the solution array
u = zeros(Nx+1, Nt+1);

% Set initial condition
u(:,1) = sin(pi*(0:Nx)*dx);
u(:,2) = u(:,1);

% Time-stepping loop
for n = 2:Nt
    for i = 2:Nx
        u(i,n+1) = 2*(1-coeff)*u(i,n) - u(i,n-1) + coeff*(u(i+1,n) + u(i-1,n));
    end
end

% Plot the solution as an animation
figure
xlabel('x')
ylabel('u(x,t)')
xlim([0, L])
ylim([-1, 1])

for n = 1:Nt
    plot(0:dx:L, u(:,n), 'LineWidth', 2);
    ylim([-1, 1]) % keep the y-limits fixed during the animation
    pause(0.01)
end

%% 
close all 
clear all

% Define spatial and temporal discretization
L = 1; % length of the domain
T = 10; % total time
Nx = 100; % number of spatial points
Nt = 1000; % number of time steps

% Define the spatial and time step sizes
dx = L/Nx;
dt = T/Nt;

% Define the wave speed and the coefficient used in the finite difference method
c = 1;
coeff = c^2*dt^2/dx^2;

% Initialize the solution array
u = zeros(Nx+1, Nt+1);

% Set initial condition
u(:,1) = sin(pi*(0:Nx)*dx);
u(:,1) = zeros(Nx+1,1);
u(:,2) = u(:,1);

% Time-stepping loop
for n = 2:Nt
    % Left Boundary condition
    u(1,n+1) = sin(2*pi*n*dt);
    for i = 2:Nx
        u(i,n+1) = 2*(1-coeff)*u(i,n) - u(i,n-1) + coeff*(u(i+1,n) + u(i-1,n));
    end
end

% Plot the solution as an animation
figure
xlabel('x')
ylabel('u(x,t)')
xlim([0, L])
ylim([-1, 1]*10)

for n = 1:2:Nt
    plot(0:dx:L, u(:,n), 'LineWidth', 2);
    ylim([-1, 1]*2 ) % keep the y-limits fixed during the animation
    pause(0.001)
end
