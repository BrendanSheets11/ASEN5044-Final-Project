clear
clc
close all

load('Rtrue.csv');
load('Qtrue.csv');
%load('orbitdeterm_finalproj_KFdata.mat');

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

%Initial values for LKF
Q = 830217.568131974*Qtrue; %1000 works good
R = 1e9*Rtrue; %1e9 works good
%Q = 0*Qtrue;

P_plus=10*[[1 0 0 0];
          [0 0.001 0 0];
          [0 0 1 0];
          [0 0 0 0.001]]; %initial state error covariance matrix
dx_plus = [0;0;0;0];

%%%Integrate non-linear EOM for true state values
% options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,x_true] = ode45(@(t,x) EOM(t,x),time,x0,options); 
num_points = 1401;
[x_true,ydata] = SimulateData(Qtrue,Rtrue,P_plus,xnom0,dt,num_points);
Xtrue = x_true(1,:);
Xdot_true = x_true(2,:);
Ytrue = x_true(3,:);
Ydot_true = x_true(4,:);
true_dx_vals = zeros([4,num_points]);
dytrue = zeros([4,num_points]);

dX_LKF = zeros([4,num_points]);
X_LKF = zeros([4,num_points]);
sigma_LKF = zeros([4,num_points]);
dyest = zeros([4,num_points]);


%%%Use discrete time linearized dynamics to estimate state values over time

for k = 0:1400

    t = dt*k; %[s] current time

    %Calculate nominal orbit state values at current time
    Xnom = r0*cos(omega0*t);
    Ynom = r0*sin(omega0*t);
    Xnom_dot = -omega0*r0*sin(omega0*t);
    Ynom_dot = omega0*r0*cos(omega0*t);

    xnom = [Xnom;Xnom_dot;Ynom;Ynom_dot]; %nominal state vector at current time
    true_dx_vals(:,k+1) = x_true(:,k+1)-xnom; %add state estimate to array

    %store values of LKF stuff
    dX_LKF(:,k+1) = dx_plus;
    X_LKF(:,k+1) = dx_plus + xnom;
    sigma_LKF(:,k+1) = 2*sqrt(diag(P_plus));
    
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
    
    %CT Perturbation Disturbance Affect matrix
    Gamma_bar =  [[0 1 0 0]; [0 0 0 1]]';

    %%%Calculate DT Perturbation Dynamics Matrices
    Ft = eye(4)+dt*A_bar; %Ahat(1:4,1:4); %DT state transition matrix at time t
    Omegat = dt*Gamma_bar;
   
    %%%Calculate Linearized Measurements
    Yk = ydata{k+2};
    
    if ~isempty(Yk)

        indices = Yk(4,:);
        H = zeros([3*length(indices),4]); %linearized measurement matrix
        Y_vect = zeros([3*length(indices),1]);
        Rk = zeros([3*length(indices),3*length(indices)]);
        t_next = t+dt;

        %Calculate nominal orbit state values at current time
        Xnom = r0*cos(omega0*t_next);
        Ynom = r0*sin(omega0*t_next);
        Xnom_dot = -omega0*r0*sin(omega0*t_next);
        Ynom_dot = omega0*r0*cos(omega0*t_next);


        %Iterate Over Each Ground Station
        for i = 1:length(indices)
    
            j = indices(i);

            ynom_i = h(j,xnom,t_next); %nominal measurement;
    
            dyi = Yk(1:3,i)-ynom_i;

            if i == 1
                dytrue(:,k+1) = [dyi;j];
            end
            
            Y_vect(3*i-2:3*i) = dyi;
            
            theta_i_t = omega_E*t_next + (j-1)*pi/6; %[rad] angle of ground station i at time t
        
            %calculate position and velocity of ground station i
            Xis = RE*cos(theta_i_t);
            Yis = RE*sin(theta_i_t);
            Xis_dot = -omega_E*RE*sin(theta_i_t);
            Yis_dot = omega_E*RE*cos(theta_i_t);
        
       
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
            H(3*i-2:3*i,:) = C_bar_i;
            Rk(3*i-2:3*i,3*i-2:3*i) = R;
            
        end
        

        %%%Time Update Step
        dx_minus = Ft*dx_plus;
        P_minus = Ft*P_plus*Ft' + Omegat*Q*Omegat';
        K = P_minus*H'/(H*P_minus*H' + Rk);

    
        %%%Measurement Update Step
        if ~isempty(indices)
            dyest(:,k+1) = [H(1:3,1:4)*dx_minus;indices(1)];
        end
        dx_plus = dx_minus + K*(Y_vect-H*dx_minus);
        P_plus = (eye(4) - K*H)*P_minus;

    else
        %%%Time Update Step
        dx_plus = Ft*dx_plus;
        P_plus = Ft*P_plus*Ft' + Omegat*Q*Omegat';
    end

end

%true_dx_vals(:,1)
% disp(true_dx_vals(:,end))
% disp(dX_LKF(:,end))
% e = x_true(:,end)-X_LKF(:,end)
% P_plus\e
% epsilon_end = e'*(P_plus\e)

%% NEES/NIS Testing
load('Rtrue.csv');
load('Qtrue.csv');

R = 1e9*Rtrue;
P_plus = 10*[[1 0 0 0];
            [0 0.001 0 0];
            [0 0 1 0];
            [0 0 0 0.001]]; %initial state error covariance matrix

alpha=0.05;
tuner = 1000;%logspace(-5,5,20);


Eps_x = zeros(length(x_true)-1,1);
Eps_bar_x = zeros(length(x_true)-1,1);

[~, ~, ~,Pk_LKF] = LKFfunc(ydata,Q,R,P_plus,dx_plus);
e_x = x_true - X_LKF;

for i=1:length(e_x)-1
    P_k = Pk_LKF(:,4*k-3:4*k);
    %Eps_x(i) = e_x(:,i)'*((diag(sigma_LKF(:,k)).^2)\e_x(:,i));
    Eps_x(i) = e_x(:,i)'*((P_k)\e_x(:,i));
    Eps_bar_x(i) = Eps_x(i)+Eps_bar_x(i);
end


r1 = chi2inv(alpha/2,4);
r2 = chi2inv(1-alpha/2,4);

num_passes = 0;

for Eps_barx_k = Eps_bar_x'
    if (Eps_barx_k >= r1) && (Eps_barx_k <= r2)
        num_passes = num_passes+1;
    end
end

average_epsilon = sum(Eps_bar_x)/length(Eps_bar_x)
proportion_failed = 1 - num_passes/num_points

if proportion_failed < alpha
    disp("KF passes chi-square test for Q:")
    disp(Q)
end

figure(3)
hold on
plot(Eps_bar_x,"bo")
plot(r1*ones(size(Eps_bar_x)),"b--")
plot(r2*ones(size(Eps_bar_x)),"b--")

%     if (N*Eps_bar_x_avg >= r1*N) && (N*Eps_bar_x_avg <= r2*N)
%         disp("KF passes chi-square test for Q:")
%         disp(Q)
%     end


% Average perturbation errors between perturbation est and LKF
XLKF_error = sum(abs(Xtrue-X_LKF(1,:)))/length(Xtrue)
YLKF_error = sum(abs(Ytrue-X_LKF(3,:)))/length(Ytrue)
XdotLKF_error = sum(abs(Xdot_true-X_LKF(2,:)))/length(Xdot_true)
YdotLKF_error = sum(abs(Ydot_true-X_LKF(4,:)))/length(Ydot_true)

figure(1)
subplot(4,1,1)
hold on
title("State vs. Time, Full Nonlinear Dynamics Simulation")
xlabel("Time (secs)")
ylabel("X (km)")
plot(time,Xtrue);
plot(time,X_LKF(1,:));
plot(time,X_LKF(1,:)+sigma_LKF(1,:),"r--");
plot(time,X_LKF(1,:)-sigma_LKF(1,:),"r--");
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("Xdot (km/s)")
plot(time,Xdot_true);
plot(time,X_LKF(2,:));
plot(time,X_LKF(2,:)+sigma_LKF(2,:),"r--");
plot(time,X_LKF(2,:)-sigma_LKF(2,:),"r--");
subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("Y (km)")
plot(time,Ytrue);
plot(time,X_LKF(3,:));
plot(time,X_LKF(3,:)+sigma_LKF(3,:),"r--");
plot(time,X_LKF(3,:)-sigma_LKF(3,:),"r--");
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("Ydot (km/s)")
plot(time,Ydot_true);
plot(time,X_LKF(4,:));
plot(time,X_LKF(4,:)+sigma_LKF(4,:),"r--");
plot(time,X_LKF(4,:)-sigma_LKF(4,:),"r--");

%{
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

%}
start = 1;
figure(4)
subplot(4,1,1)
hold on
title("Approximate Perturbation Errors vs. Time")
xlabel("Time (secs)")
ylabel("\deltaX (km)")
plot(time(start:end),true_dx_vals(1,start:end),"b");
plot(time(start:end),dX_LKF(1,start:end),"r");
plot(time(start:end),dX_LKF(1,start:end)+sigma_LKF(1,start:end),"r--");
plot(time(start:end),dX_LKF(1,start:end)-sigma_LKF(1,start:end),"r--");
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("\deltaXdot (km/s)")

plot(time(start:end),true_dx_vals(2,start:end),"b");
plot(time(start:end),dX_LKF(2,start:end),"r");
plot(time(start:end),dX_LKF(2,start:end)+sigma_LKF(2,start:end),"r--");
plot(time(start:end),dX_LKF(2,start:end)-sigma_LKF(2,start:end),"r--");

subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("\deltaY (km)")
plot(time(start:end),true_dx_vals(3,start:end),"b");
plot(time(start:end),dX_LKF(3,start:end),"r");
plot(time(start:end),dX_LKF(3,start:end)+sigma_LKF(3,start:end),"r--");
plot(time(start:end),dX_LKF(3,start:end)-sigma_LKF(3,start:end),"r--");
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("\deltaYdot (km/s)")
plot(time(start:end),true_dx_vals(4,start:end),"b");
plot(time(start:end),dX_LKF(4,start:end),"r");
plot(time(start:end),dX_LKF(4,start:end)+sigma_LKF(4,start:end),"r--");
plot(time(start:end),dX_LKF(4,start:end)-sigma_LKF(4,start:end),"r--");

% start = 100;
% 
figure(5)
subplot(4,1,1)
hold on
title("Approximate Perturbation Errors vs. Time")
xlabel("Time (secs)")
ylabel("\deltaX (km)")
plot(time(start:end),true_dx_vals(1,start:end)-dX_LKF(1,start:end),"b");
plot(time(start:end),sigma_LKF(1,start:end),"r--");
plot(time(start:end),-sigma_LKF(1,start:end),"r--");
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("\deltaXdot (km/s)")

plot(time(start:end),true_dx_vals(2,start:end)-dX_LKF(2,start:end),"b");
plot(time(start:end),sigma_LKF(2,start:end),"r--");
plot(time(start:end),-sigma_LKF(2,start:end),"r--");

subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("\deltaY (km)")
plot(time(start:end),true_dx_vals(3,start:end)-dX_LKF(3,start:end),"b");
plot(time(start:end),sigma_LKF(3,start:end),"r--");
plot(time(start:end),-sigma_LKF(3,start:end),"r--");
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("\deltaYdot (km/s)")

plot(time(start:end),true_dx_vals(4,start:end)-dX_LKF(4,start:end),"b");
plot(time(start:end),sigma_LKF(4,start:end),"r--");
plot(time(start:end),-sigma_LKF(4,start:end),"r--");

%Plot true measurements
figure(2)
subplot(4,1,1)
hold on
scatter(time,dyest(1,:))
scatter(time,dytrue(1,:))
title("Full Nonlinear Model Data Simulation")
xlabel("Time (secs)")
ylabel("\rho^i (km)")
subplot(4,1,2)
scatter(time,dyest(2,:))
hold on
scatter(time,dytrue(2,:))
xlabel("Time (secs)")
ylabel("$\dot{\rho}^i$ (km/s)",Interpreter="latex")
subplot(4,1,3)
scatter(time,dyest(3,:))
hold on
scatter(time,dytrue(3,:))
xlabel("Time (secs)")
ylabel("\phi^i (rads)")
subplot(4,1,4)
scatter(time,dyest(4,:))
hold on
scatter(time,dytrue(4,:))
xlabel("Time (secs)")
ylabel("Visible Station ID")


%{
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

%}

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

function valid_meas = all_measurements(x,t)
    
    valid_meas = [];

    for i = 1:12
        
        omega_E = 2*pi/86400; %[rad/s] angular velocity of the earth
        theta_i_t = omega_E*t + (i-1)*pi/6; %[rad] angle of ground station i at time t

        meas_i = h(i,x,t);

        theta_i = atan2(sin(theta_i_t),cos(theta_i_t));

        angle_diff = theta_i - meas_i(3);
        angle_diff = mod(angle_diff + pi,2*pi) - pi; %[rad] wrapped angle difference between two angles 
        
        %store observations if angle difference is less than 90 degrees
        if abs(angle_diff) < pi/2
            valid_meas = [valid_meas, [meas_i;i]];
        end
    end

end

function [X_sim,Y_sim] = SimulateData(Qtrue,Rtrue,P0,Xnom0,dt,num_points)
   
   cholP = chol(P0,"lower");
   noisy_state = cholP*randn([4,1]) + Xnom0;
   X_sim = zeros([4,num_points]);
   Y_sim = cell(1,num_points+1);

   cholQ = chol(Qtrue,"lower");
   cholR = chol(Rtrue,"lower");

   for i = 1:num_points
        

        %generate noisy measurements
        true_meas = all_measurements(noisy_state,dt*(i-1));
        
        if ~isempty(true_meas)
            noisy_meas = zeros(size(true_meas));
            for j = 1:length(true_meas(1,:))
                rand_j = randn([3,1]);
                measurement_noise_j = cholR*rand_j;
                noisy_meas(:,j) = true_meas(:,j) + [measurement_noise_j;0];
            end
        else
            noisy_meas = [];
        end
        
        %store previous state and measurement
        X_sim(:,i) = noisy_state; 
        Y_sim{i} = noisy_meas;

        tspan = [0 dt];

        %%%Integrate non-linear EOM for true state values
        rand1 = randn([2,1]);
        process_noise = cholQ*rand1;
        options = odeset('RelTol',1e-8,'AbsTol',1e-8);
        [~,trajectory] = ode45(@(t,x) noisy_EOM(t,x,process_noise),tspan,noisy_state,options); 
        
        %save next state
        noisy_state = trajectory(end,:);
   end

   %generate final state measurement
   true_meas = all_measurements(noisy_state,dt*(num_points+1));
    
   if ~isempty(true_meas)
       noisy_meas = zeros(size(true_meas));
       for j = 1:length(true_meas(1,:))
           rand_j = randn([3,1]);
           measurement_noise_j = cholR*rand_j;
           noisy_meas(:,j) = true_meas(:,j) + [measurement_noise_j;0];
       end
   else
       noisy_meas = [];
   end
    
    %store previous state and measurement
    Y_sim{num_points+1} = noisy_meas;
end

function Xdot = noisy_EOM(t,X,w)
    mu = 398600; %[km^3/s^2] gravitational parameter
    r = sqrt(X(1)^2+X(3)^2);
    Xdot = [X(2); -mu*X(1)/r^3 + w(1); X(4); -mu*X(3)/r^3 + w(2)];
end
