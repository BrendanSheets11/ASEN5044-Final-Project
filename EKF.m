function  [x_true, ydata, X_EKF, sigma_EKF,Pk_EKF,eps_y] = EKF(Q,R,P0)

load('Rtrue.csv');
load('Qtrue.csv');

%Define important constants
r0 = 6678; %[km] nominal orbit radius
mu = 398600; %[km^3/s^2] gravitational parameter
omega0 = sqrt(mu/r0^3); %[rad/s] nominal orbit velocity
dt = 10; %[s] simulation time step
time = 0:dt:14000; %[s] simulation time step array
RE = 6378; %[km] radius of the earth
omega_E = 2*pi/86400; %[rad/s] angular velocity of the earth

dx = [0;0.075;0;-0.021]; %initial state perturbation
Xnom0 = [r0;0;0;omega0*r0]; %initial nominal state




%Initial values for LKF

P_plus = P0;
x_plus =  Xnom0; %initial state

%%%Integrate non-linear EOM for true state values
% options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,x_true] = ode45(@(t,x) EOM(t,x),time,x0,options); 
num_points = 1401;
[x_true,ydata] = SimulateData(Qtrue,Rtrue,P_plus,Xnom0,dt,num_points);
Xtrue = x_true(1,:);
Xdot_true = x_true(2,:);
Ytrue = x_true(3,:);
Ydot_true = x_true(4,:);
ytrue = zeros([4,num_points]);

X_EKF = zeros([4,num_points]);
sigma_EKF = zeros([4,num_points]);
yest = zeros([4,num_points]);
eps_y = zeros([1,length(ydata)]);
Pk_EKF = zeros([4,4*num_points]);

%%%Use discrete time linearized dynamics to estimate state values over time

for k = 0:1400

    t = dt*k; %[s] current time

    %store values of LKF stuff
    X_EKF(:,k+1) = x_plus;
    sigma_EKF(:,k+1) = 2*sqrt(diag(P_plus));
    Pk_EKF(:,4*k+1:4*k+4) = P_plus;
    %Calculate partial derivatives for dynamics jacobian evaluated on
    %nominal trajectory at the current time
    Xnow = x_plus(1);
    Xnow_dot = x_plus(2);
    Ynow = x_plus(3);
    Ynow_dot = x_plus(4);
    rnow = sqrt(Xnow^2+Ynow^2);
    dxdot_dx = mu*(2*Xnow^2 - Ynow^2)/(rnow^5);
    dxdot_dy = 3*mu*Xnow*Ynow/(rnow^5);
    dydot_dy = mu*(2*Ynow^2 - Xnow^2)/(rnow^5);
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


        %Iterate Over Each Ground Station
        for i = 1:length(indices)
    
            j = indices(i);
    
            yi = Yk(1:3,i);

            if i == 1
                ytrue(:,k+1) = [yi;j];
            end
            
            Y_vect(3*i-2:3*i) = yi;
            
            theta_i_t = omega_E*t_next + (j-1)*pi/6; %[rad] angle of ground station i at time t
        
            %calculate position and velocity of ground station i
            Xis = RE*cos(theta_i_t);
            Yis = RE*sin(theta_i_t);
            Xis_dot = -omega_E*RE*sin(theta_i_t);
            Yis_dot = omega_E*RE*cos(theta_i_t);
        
       
            %Calculate Linearized Observation Matrix
            rho = sqrt((Xnow-Xis)^2+(Ynow-Yis)^2);
    
            drho_dx = (Xnow-Xis)/rho;
            drho_dy = (Ynow-Yis)/rho;
            
            drhodot_dx = (Ynow-Yis)*((Xnow_dot-Xis_dot)*(Ynow-Yis) - (Xnow-Xis)*(Ynow_dot-Yis_dot))/rho^3;
            drhodot_dy = (Xnow-Xis)*((Xnow-Xis)*(Ynow_dot-Yis_dot) - (Xnow_dot-Xis_dot)*(Ynow-Yis))/rho^3;
            drhodot_dxdot = (Xnow-Xis)/rho;
            drhodot_dydot = (Ynow-Yis)/rho;
    
            dphi_dx = (Yis-Ynow)/rho^2;
            dphi_dy = (Xnow-Xis)/rho^2;
            
            C_bar_i = [[drho_dx    0             drho_dy    0];
                       [drhodot_dx drhodot_dxdot drhodot_dy drhodot_dydot];
                       [dphi_dx    0             dphi_dy    0]];
            H(3*i-2:3*i,:) = C_bar_i;
            Rk(3*i-2:3*i,3*i-2:3*i) = R;
            
        end
        

        %%%Time Update Step
        %dx_minus = Ft*dx_plus;
        options = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [~,x_traj] = ode45(@(t,x) noisy_EOM(t,x,[0;0]),[0 dt],x_plus,options); 
        x_minus = x_traj(end,:)';

        P_minus = Ft*P_plus*Ft' + Omegat*Q*Omegat';
        K = P_minus*H'/(H*P_minus*H' + Rk);
      
    
        %%%Measurement Update Step
        Yhat = zeros([3*length(indices),1]); %estimated measurements of current state

        %calculate what the nonlinear measurements of the current state
        %estimate should be
        if ~isempty(indices)
             
            %Iterate Over Each Ground Station
            for i = 1:length(indices)
        
                j = indices(i);

                yhat_i = h(j,x_plus,t_next); %nominal measurement;
                Yhat(3*i-2:3*i) = yhat_i;
            end
        end
        
        x_plus = x_minus + K*(Y_vect-Yhat);
        P_plus = (eye(4) - K*H)*P_minus;

        e_y = Y_vect-Yhat;
        eps_y(k+1) = e_y'*((H*P_minus*H' + Rk/1e6)\e_y);
    else
        %%%Time Update Step
        %dx_plus = Ft*dx_plus;
        [~,x_traj] = ode45(@(t,x) noisy_EOM(t,x,[0;0]),[0 dt],x_plus,options); 
        x_plus = x_traj(end,:)';
        P_plus = Ft*P_plus*Ft' + Omegat*Q*Omegat';

    end

end

%{

% Average perturbation errors between perturbation est and LKF
XEKF_error = sum(abs(Xtrue-X_EKF(1,:)))/length(Xtrue)
YEKF_error = sum(abs(Ytrue-X_EKF(3,:)))/length(Ytrue)
XdotEKF_error = sum(abs(Xdot_true-X_EKF(2,:)))/length(Xdot_true)
YdotEKF_error = sum(abs(Ydot_true-X_EKF(4,:)))/length(Ydot_true)

figure(1)
subplot(4,1,1)
hold on
title("State vs. Time, Full Nonlinear Dynamics Simulation")
xlabel("Time (secs)")
ylabel("X (km)")
plot(time,Xtrue);
plot(time,X_EKF(1,:));
plot(time,X_EKF(1,:)+sigma_EKF(1,:),"r--");
plot(time,X_EKF(1,:)-sigma_EKF(1,:),"r--");
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("Xdot (km/s)")
plot(time,Xdot_true);
plot(time,X_EKF(2,:));
plot(time,X_EKF(2,:)+sigma_EKF(2,:),"r--");
plot(time,X_EKF(2,:)-sigma_EKF(2,:),"r--");
subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("Y (km)")
plot(time,Ytrue);
plot(time,X_EKF(3,:));
plot(time,X_EKF(3,:)+sigma_EKF(3,:),"r--");
plot(time,X_EKF(3,:)-sigma_EKF(3,:),"r--");
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("Ydot (km/s)")
plot(time,Ydot_true);
plot(time,X_EKF(4,:));
plot(time,X_EKF(4,:)+sigma_EKF(4,:),"r--");
plot(time,X_EKF(4,:)-sigma_EKF(4,:),"r--");

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

figure(5)
subplot(4,1,1)
hold on
title("Approximate Errors vs. Time")
xlabel("Time (secs)")
ylabel("\deltaX (km)")
plot(time(start:end),x_true(1,start:end)-X_EKF(1,start:end),"b");
% plot(time(start:end),sigma_EKF(1,start:end),"r--");
% plot(time(start:end),-sigma_EKF(1,start:end),"r--");
subplot(4,1,2)
hold on
xlabel("Time (secs)")
ylabel("\deltaXdot (km/s)")

plot(time(start:end),x_true(2,start:end)-X_EKF(2,start:end),"b");
% plot(time(start:end),sigma_EKF(2,start:end),"r--");
% plot(time(start:end),-sigma_EKF(2,start:end),"r--");

subplot(4,1,3)
hold on
xlabel("Time (secs)")
ylabel("\deltaY (km)")
plot(time(start:end),x_true(3,start:end)-X_EKF(3,start:end),"b");
% plot(time(start:end),sigma_EKF(3,start:end),"r--");
% plot(time(start:end),-sigma_EKF(3,start:end),"r--");
subplot(4,1,4)
hold on
xlabel("Time (secs)")
ylabel("\deltaYdot (km/s)")

plot(time(start:end),x_true(4,start:end)-X_EKF(4,start:end),"b");
% plot(time(start:end),sigma_EKF(4,start:end),"r--");
% plot(time(start:end),-sigma_EKF(4,start:end),"r--");

%Plot true measurements
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
%}
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

function [X_sim,Y_sim] = SimulateData(Qtrue,Rtrue,P0,Xnow0,dt,num_points)
   
   cholP = chol(P0,"lower");
   noisy_state = cholP*randn([4,1]) + Xnow0;
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
