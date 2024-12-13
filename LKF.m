clear
%clc
close all

load('Rtrue.csv');
load('Qtrue.csv');
load('orbitdeterm_finalproj_KFdata.mat');

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
Qtuning = 1;
Ptuning = 1e-1;
% Q = zeros([4,4]);
% Q(2,2) = Qtrue(1,1); %dynamics noise covariance matrix
% Q(4,4) = Qtrue(2,2);
Q = Qtuning*Qtrue;
%P_plus = 1e6*eye(4);
P_plus = Ptuning*[[10 0 0 0];
                 [0 1 0 0];
                 [0 0 10 0];
                 [0 0 0 1]]; %initial state error covariance matrix
dx_plus = dx;

%%%Integrate non-linear EOM for true state values
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x_true] = ode45(@(t,x) EOM(t,x),time,x0,options); 
Xtrue = x_true(:,1);
Xdot_true = x_true(:,2);
Ytrue = x_true(:,3);
Ydot_true = x_true(:,4);



dyest = zeros([4,length(ydata)]);
dytrue = zeros([4,length(ydata)]);

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
true_dx_vals = zeros([4,1401]);
meas2 = zeros([12,3,length(Xtrue)]); %stores list of all valid linearized measurements for all ground stations over all time



%% NEES/NIS Testing

alpha=0.05;
tuner = logspace(-5,5,20);
for j=1:length(tuner)
    
    Q = tuner(j);

    Eps_x = zeros(length(x_true)-1,1);
    Eps_bar_x = zeros(length(x_true)-1,1);
    
    for N=1:30
        [dX_LKF, X_LKF, sigma_LKF,P_k] = LKFfunc(ydata,Q,Rtrue,P_plus,dx);
        %e_x = X_sim - X_LKF
        e_x = x_true' - X_LKF;
        for i=1:length(e_x)-1
            Eps_x(i) = e_x(:,i)'*((P_k)\e_x(:,i));
            Eps_bar_x(i) = Eps_x(i)+Eps_bar_x(i);
        end
        
    end
    Eps_bar_x = Eps_bar_x/N;
    Eps_bar_x_avg = sum(Eps_bar_x)/length(Eps_bar_x);
    
    Eps_Q_test(j) = Eps_bar_x_avg;

    r1 = chi2inv(alpha/2,N*4)./N;
    r2 = chi2inv(1-alpha/2,N*4)./N;
    
    if (N*Eps_bar_x_avg >= r1*N) && (N*Eps_bar_x_avg <= r2*N)
        disp("KF passes chi-square test for Q:")
        disp(Q)
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
Xpos_error = sum(abs(Xtrue-Xest'))/length(Xtrue);
Ypos_error = sum(abs(Ytrue-Yest'))/length(Ytrue);
Xvel_error = sum(abs(Xdot_true-Xdot_est'))/length(Xdot_true);
Yvel_error = sum(abs(Ydot_true-Ydot_est'))/length(Ydot_true);

% Average perturbation errors between perturbation est and LKF
% XLKF_error = sum(abs(Xtrue-X_LKF(1,:)'))/length(Xtrue)
% YLKF_error = sum(abs(Ytrue-X_LKF(3,:)'))/length(Ytrue)
% XdotLKF_error = sum(abs(Xdot_true-X_LKF(2,:)'))/length(Xdot_true)
% YdotLKF_error = sum(abs(Ydot_true-X_LKF(4,:)'))/length(Ydot_true)



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

% figure(4)
% subplot(4,1,1)
% hold on
% title("Linearized Approximate Perturbations vs. Time")
% xlabel("Time (secs)")
% ylabel("\deltaX (km)")
% plot(time,true_dx_vals(1,:));
% plot(time,dX_LKF(1,:));
% plot(time,dX_LKF(1,:)+sigma_LKF(1,:));
% plot(time,dX_LKF(1,:)-sigma_LKF(1,:));
% subplot(4,1,2)
% hold on
% xlabel("Time (secs)")
% ylabel("\deltaXdot (km/s)")
% start = 100;
% % plot(time(start:end),true_dx_vals(2,start:end));
% % plot(time(start:end),dX_LKF(2,start:end));
% % plot(time(start:end),dX_LKF(2,start:end)+sigma_LKF(2,start:end));
% % plot(time(start:end),dX_LKF(2,start:end)-sigma_LKF(2,start:end));
% 
% plot(time(start:end),true_dx_vals(2,start:end)-dX_LKF(2,start:end),"b");
% plot(time(start:end),sigma_LKF(2,start:end),"r--");
% plot(time(start:end),-sigma_LKF(2,start:end),"r--");
% 
% subplot(4,1,3)
% hold on
% xlabel("Time (secs)")
% ylabel("\deltaY (km)")
% plot(time,true_dx_vals(3,:));
% plot(time,dX_LKF(3,:));
% plot(time,dX_LKF(3,:)+sigma_LKF(3,:));
% plot(time,dX_LKF(3,:)-sigma_LKF(3,:));
% subplot(4,1,4)
% hold on
% xlabel("Time (secs)")
% ylabel("\deltaYdot (km/s)")
% 
% plot(time(start:end),true_dx_vals(4,start:end));
% plot(time(start:end),dX_LKF(4,start:end));
% plot(time(start:end),dX_LKF(4,start:end)+sigma_LKF(4,start:end));
% plot(time(start:end),dX_LKF(4,start:end)-sigma_LKF(4,start:end));

start = 100;

figure(4)
subplot(4,1,1)
hold on
title("Linearized Approximate Perturbations vs. Time")
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

% figure()
% %Plot true measurements
% figure(2)
% subplot(4,1,1)
% hold on
% scatter(time,dyest(1,:))
% scatter(time,dytrue(1,:))
% title("Full Nonlinear Model Data Simulation")
% xlabel("Time (secs)")
% ylabel("\rho^i (km)")
% subplot(4,1,2)
% scatter(time,dyest(2,:))
% hold on
% scatter(time,dytrue(2,:))
% xlabel("Time (secs)")
% ylabel("$\dot{\rho}^i$ (km/s)",Interpreter="latex")
% subplot(4,1,3)
% scatter(time,dyest(3,:))
% hold on
% scatter(time,dytrue(3,:))
% xlabel("Time (secs)")
% ylabel("\phi^i (rads)")
% subplot(4,1,4)
% scatter(time,dyest(4,:))
% hold on
% scatter(time,dytrue(4,:))
% xlabel("Time (secs)")
% ylabel("Visible Station ID")


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
