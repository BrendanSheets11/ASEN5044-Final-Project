close all
clc
%% NEES/NIS Testing
load('Rtrue.csv');
load('Qtrue.csv');

R = 1e9*Rtrue;
P_plus = 10*[[1 0 0 0];
            [0 0.001 0 0];
            [0 0 1 0];
            [0 0 0 0.001]]; %initial state error covariance matrix

alpha=0.05;
tuner = logspace(6.14,6.145,10);
for j=1:length(tuner)
    
    Q = tuner(j)*Qtrue;

    Eps_x = zeros(length(x_true)-1,1);
    Eps_bar_x = zeros(length(x_true)-1,1);
    
    for N=1:30
        [x_true,ydata] = SimulateData(Qtrue,Rtrue,P_plus,xnom0,dt,num_points);
        [dX_LKF, X_LKF, sigma_LKF,Pk_LKF] = LKFfunc(ydata,Q,R,P_plus,dx_plus);
        e_x = x_true - X_LKF;

        for i=1:length(e_x)-1
            P_k = Pk_LKF(:,4*k-3:4*k);
            Eps_x(i) = e_x(:,i)'*((P_k)\e_x(:,i));
            Eps_bar_x(i) = Eps_x(i)+Eps_bar_x(i);
        end
        
    end
    Eps_bar_x = Eps_bar_x/N;
    Eps_bar_x_avg = sum(Eps_bar_x)/length(Eps_bar_x)
    
    Eps_Q_test(j) = Eps_bar_x_avg;

    r1 = chi2inv(alpha/2,N*4)./N;
    r2 = chi2inv(1-alpha/2,N*4)./N;

    num_passes = 0;

    for Eps_barx_k = Eps_bar_x'
        if (N*Eps_barx_k >= r1*N) && (N*Eps_barx_k <= r2*N)
            num_passes = num_passes+1;
        end
    end

    proportion_failed = 1 - num_passes/num_points;

    if proportion_failed < alpha
        disp("KF passes chi-square test for Q:")
        disp(Q)
    end
    
    figure(j)
    hold on
    plot(Eps_bar_x,"bo")
    plot(r1*ones(size(Eps_bar_x)),"r--")
    plot(r2*ones(size(Eps_bar_x)),"r--")
    ylim([0 100]);
    hold off
    
%     if (N*Eps_bar_x_avg >= r1*N) && (N*Eps_bar_x_avg <= r2*N)
%         disp("KF passes chi-square test for Q:")
%         disp(Q)
%     end
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

function [X_sim,Y_sim] = SimulateData(Qtrue,Rtrue,P0,Xnom0,dt,num_points)
   
   cholP = chol(P0,"lower");
   noisy_state = cholP*randn([4,1]) + Xnom0;
   X_sim = zeros([4,num_points]);
   Y_sim = cell(1,num_points);

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
end

function Xdot = noisy_EOM(t,X,w)
    mu = 398600; %[km^3/s^2] gravitational parameter
    r = sqrt(X(1)^2+X(3)^2);
    Xdot = [X(2); -mu*X(1)/r^3 + w(1); X(4); -mu*X(3)/r^3 + w(2)];
end

