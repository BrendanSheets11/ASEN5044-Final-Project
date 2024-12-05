clear
clc
close all

%Define important constants
r0 = 6678; %[km] nominal orbit radius
mu = 398600; %[km^3/s^2] gravitational parameter
omega0 = sqrt(mu/r0^3); %[rad/s] nominal orbit velocity
dt = 10; %[s] simulation time step
time = 0:dt:14000; %[s] simulation time step array
RE = 6378; %[km] radius of the earth
omega_E = 2*pi/86400; %[rad/s] angular velocity of the earth

dx = [0;0.075;0;-0.021]; %initial state perturbation
xnom0 = [r0;0;0;omega0*r0]; %initial nominal state
x0 = xnom0+dx; %initial state

%%%Integrate non-linear EOM for true state values
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x_true] = ode45(@(t,x) EOM(t,x),time,x0,options); 
Xtrue = x_true(:,1);
Xdot_true = x_true(:,2);
Ytrue = x_true(:,3);
Ydot_true = x_true(:,4);

%calculate non-linear measurements based on non-linear true state
meas = zeros([12,3,length(Xtrue)]); %stores list of all valid measurements for all ground stations over all time
all_true_angles = zeros([12,length(Xtrue)]); %stores list of all angles for all ground stations over all time

for k = 0:length(Xtrue)-1 %iterate over all time steps

    for i = 1:12 %iterate over all ground stations

        yi = h(i,x_true(k+1,:),k*dt); %generate measurement of ground station at current time
        all_true_angles(i,k+1) = yi(3); %store true value of phi for later use

        theta_i = omega_E*k*dt + (i-1)*pi/6; %[rad] angle of ground station i at time t
        theta_i = atan2(sin(theta_i),cos(theta_i));

        angle_diff = theta_i - yi(3);
        angle_diff = mod(angle_diff + pi,2*pi) - pi; %[rad] wrapped angle difference between two angles 
        
        %store observations if angle difference is less than 90 degrees
        if abs(angle_diff) < pi/2
            meas(i,:,k+1) = yi;
        end
    end
end


%%%Use discrete time linearized dynamics to estimate state values over time
dx_vals = zeros([4,1401]);
x_vals = zeros([4,1401]);
meas2 = zeros([12,3,length(Xtrue)]); %stores list of all valid linearized measurements for all ground stations over all time

for k = 0:1400

    dx_vals(:,k+1) = dx; %add perturbation to array

    t = dt*k; %[s] current time

    %Calculate nominal orbit state values at current time
    Xnom = r0*cos(omega0*t);
    Ynom = r0*sin(omega0*t);
    Xnom_dot = -omega0*r0*sin(omega0*t);
    Ynom_dot = omega0*r0*cos(omega0*t);

    xnom = [Xnom;Xnom_dot;Ynom;Ynom_dot]; %nominal state vector at current time
    x = xnom+dx; %linearized estimated state at current time
    x_vals(:,k+1) = x; %add state estimate to array
    
    %Calculate partial derivatives for dynamics jacobian evaluated on
    %nominal trajectory at the current time
    dxdot_dx = mu*(2*Xnom^2 - Ynom^2)/(r0^5);
    dxdot_dy = 3*mu*Xnom*Ynom/(r0^5);
    dydot_dy = mu*(2*Ynom^2 - Xnom^2)/(r0^5);
    dydot_dx = dxdot_dy;
    
    %CT Perturbation Dynamics State transition matrix at current time
    A_bar = [[0       1 0        0];
            [dxdot_dx 0 dxdot_dy 0];
            [0        0 0        1];
            [dydot_dx 0 dydot_dy 0]];
    
    %CT Perturbation Control Affect matrix
    B_bar = [[0 1 0 0]; [0 0 0 1]]';
    
    %CT Perturbation Disturbance Affect matrix
    Gamma_bar = B_bar;

    %%%Calculate DT Perturbation Dynamics Matrices
    Ahat = [[A_bar, B_bar];[zeros([2,6])]];
    eAhat = expm(dt*Ahat);

    Ft = eAhat(1:4,1:4); %DT state transition matrix at time t
    Gt = eAhat(1:4,5:6); %DT control affect matrix at time t

    dx = Ft*dx; %Update the state perturbation
    
   
    %%%Calculate Linearized Measurements
    
    %Iterate Over Each Ground Station
    for i = 1:12

        ynom_i = h(i,xnom,t); %nominal measurement
        
        theta_i_t = omega_E*t + (i-1)*pi/6; %[rad] angle of ground station i at time t
    
        %calculate position and velocity of ground station i
        Xis = RE*cos(theta_i_t);
        Yis = RE*sin(theta_i_t);
        Xis_dot = -omega_E*RE*sin(theta_i_t);
        Yis_dot = omega_E*RE*cos(theta_i_t);
    
        theta_i = atan2(Yis,Xis);
        angle_diff = theta_i - ynom_i(3); %all_true_angles(i,k+1); %use true phi value to calculate if station is observable
        angle_diff = mod(angle_diff + pi,2*pi) - pi; %[rad] wrapped angle difference between two angles 
        
        %store observations if angle difference is less than 90 degrees
        if abs(angle_diff) < pi/2
    
            %Calculate Linearized Observation Matrix
            rho = sqrt((Xnom-Xis)^2+(Ynom-Yis)^2);
    
            drho_dx = (Xnom-Xis)/rho;
            drho_dy = (Ynom-Yis)/rho;
            
            drhodot_dx = (Ynom-Yis)*((Xnom_dot-Xis_dot)*(Ynom-Yis) - (Xnom-Xis)*(Ynom_dot-Yis_dot))/rho^3;
            drhodot_dy = (Xnom-Xis)*((Xnom-Xis)*(Ynom_dot-Yis_dot) - (Xnom_dot-Xis_dot)*(Ynom-Yis))/rho^3;
            drhodot_dxdot = (Xnom-Xis)/rho;
            drhodot_dydot = (Ynom-Yis)/rho;
    
            dphi_dx = (Yis-Ynom)/rho^2;
            dphi_dy = (Xnom-Xis)/rho^2;
            
            C_bar_i = [[drho_dx    0             drho_dy    0];
                       [drhodot_dx drhodot_dxdot drhodot_dy drhodot_dydot];
                       [dphi_dx    0             dphi_dy    0]];
    
            dy_i = C_bar_i*dx; %measurement perturbation

            y_i = ynom_i + dy_i; %full linearized measurement prediction
            y_i(3) = mod(y_i(3)+pi,2*pi)-pi;
            meas2(i,:,k+1) = y_i;
        end
    end
end

%extract estimated state values
Xest = x_vals(1,:);
Xdot_est = x_vals(2,:);
Yest = x_vals(3,:);
Ydot_est = x_vals(4,:);

%extract estimated state perturbations
dXest = dx_vals(1,:);
dXdot_est = dx_vals(2,:);
dYest = dx_vals(3,:);
dYdot_est = dx_vals(4,:);


% Finding average state errors
Xpos_error = sum(abs(Xtrue-Xest'))/length(Xtrue)
Ypos_error = sum(abs(Ytrue-Yest'))/length(Ytrue)
Xvel_error = sum(abs(Xdot_true-Xdot_est'))/length(Xdot_true)
Yvel_error = sum(abs(Ydot_true-Ydot_est'))/length(Ydot_true)


figure(1)
subplot(4,1,1)
hold on
title("State vs. Time, Full Nonlinear Dynamics Simulation")
xlabel("Time (secs)")
ylabel("X (km)")
plot(time,Xtrue);
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("Xdot (km/s)")
plot(time,Xdot_true);
subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("Y (km)")
plot(time,Ytrue);
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("Ydot (km/s)")
plot(time,Ydot_true);

%Plot true measurements
figure(2)
subplot(4,1,1)
title("Full Nonlinear Model Data Simulation")
hold on
xlabel("Time (secs)")
ylabel("\rho^i (km)")
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("$\dot{\rho}^i$ (km/s)",Interpreter="latex")
subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("\phi^i (rads)")
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("Visible Station ID")
for i = 1:12
    rho_i = meas(i,1,:);
    rhodot_i = meas(i,2,:);
    phi_i = meas(i,3,:);
    t_i = time(rho_i~=0);
    rhodot_i = rhodot_i(rho_i ~= 0);
    phi_i = phi_i(rho_i ~= 0);
    rho_i = rho_i(rho_i ~= 0);
    new_rho = zeros([1,length(rho_i)]);
    new_rhodot = zeros([1,length(rhodot_i)]);
    new_phi = zeros([1,length(phi_i)]);
    for j = 1:length(rho_i)
        new_rho(j) = rho_i(j);
        new_rhodot(j) = rhodot_i(j);
        new_phi(j) = phi_i(j);
    end
    subplot(4,1,1)
    hold on
    scatter(t_i,new_rho,"o")
    subplot(4,1,2)
    hold on
    scatter(t_i,new_rhodot,"o")
    subplot(4,1,3)
    hold on
    scatter(t_i,new_phi,"o")
    subplot(4,1,4)
    hold on
    scatter(t_i,i*ones(length(t_i)),"^")
end

figure(3)
subplot(4,1,1)
hold on
title("States vs Time, Linearized Approximate Dynamics Simulation")
xlabel("Time (secs)")
ylabel("X (km)")
plot(time,Xest);
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("Xdot (km/s)")
plot(time,Xdot_est);
subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("Y (km)")
plot(time,Yest);
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("Ydot (km/s)")
plot(time,Ydot_est);

figure(4)
subplot(4,1,1)
hold on
title("Linearized Approximate Perturbations vs. Time")
xlabel("Time (secs)")
ylabel("\deltaX (km)")
plot(time,dXest);
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("\deltaXdot (km/s)")
plot(time,dXdot_est);
subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("\deltaY (km)")
plot(time,dYest);
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("\deltaYdot (km/s)")
plot(time,dYdot_est);

%Plot estimated measurements
figure(5)
subplot(4,1,1)
title("Approximate Linearized Model Data Simulation")
hold on
xlabel("Time (secs)")
ylabel("\rho^i (km)")
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("$\dot{\rho}^i$ (km/s)",Interpreter="latex")
subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("\phi^i (rads)")
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("Visible Station ID")
for i = 1:12
    rho_i = meas2(i,1,:);
    rhodot_i = meas2(i,2,:);
    phi_i = meas2(i,3,:);
    t_i = time(rho_i~=0);
    rhodot_i = rhodot_i(rho_i ~= 0);
    phi_i = phi_i(rho_i ~= 0);
    rho_i = rho_i(rho_i ~= 0);
    new_rho = zeros([1,length(rho_i)]);
    new_rhodot = zeros([1,length(rhodot_i)]);
    new_phi = zeros([1,length(phi_i)]);
    for j = 1:length(rho_i)
        new_rho(j) = rho_i(j);
        new_rhodot(j) = rhodot_i(j);
        new_phi(j) = phi_i(j);
    end
    subplot(4,1,1)
    hold on
    scatter(t_i,new_rho,"o")
    subplot(4,1,2)
    hold on
    scatter(t_i,new_rhodot,"o")
    subplot(4,1,3)
    hold on
    scatter(t_i,new_phi,"o")
    subplot(4,1,4)
    hold on
    scatter(t_i,i*ones(length(t_i)),"^")

end

function Xdot = EOM(t,X)
    mu = 398600; %[km^3/s^2] gravitational parameter
    r = sqrt(X(1)^2+X(3)^2);
    Xdot = [X(2); -mu*X(1)/r^3; X(4); -mu*X(3)/r^3];
end

function meas = h(i,x,t)
    %extract state info
    X = x(1);
    Xdot = x(2);
    Y = x(3);
    Ydot = x(4);

    %calculate position and velocity of ground station i
    RE = 6378; %[km] radius of the earth
    omega_E = 2*pi/86400; %[rad/s] angular velocity of the earth
    theta_i_t = omega_E*t + (i-1)*pi/6; %[rad] angle of ground station i at time t
    Xis = RE*cos(theta_i_t);
    Yis = RE*sin(theta_i_t);
    Xis_dot = -omega_E*RE*sin(theta_i_t);
    Yis_dot = omega_E*RE*cos(theta_i_t);
    
    %calculate measurements
    rho = sqrt((X-Xis)^2+(Y-Yis)^2);
    rho_dot = ( (X-Xis)*(Xdot-Xis_dot) + (Y-Yis)*(Ydot-Yis_dot))/rho;
    phi = atan2(Y-Yis,X-Xis);
    meas = [rho;rho_dot;phi];
end
