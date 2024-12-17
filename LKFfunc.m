function [dX_LKF, X_LKF, sigma_LKF,Pk_LKF,eps_y] = LKFfunc(ydata,Q,Rtrue,P0,dx0)
r0 = 6678; %[km] nominal orbit radius
mu = 398600; %[km^3/s^2] gravitational parameter
omega0 = sqrt(mu/r0^3); %[rad/s] nominal orbit velocity
dt = 10; %[s] simulation time step
RE = 6378; %[km] radius of the earth
omega_E = 2*pi/86400; %[rad/s] angular velocity of the earth
%%% Not sure about R

dX_LKF = zeros([4,length(ydata)]);
X_LKF = zeros([4,length(ydata)]);
sigma_LKF = zeros([4,length(ydata)]);
Pk_LKF = zeros([4,4*length(ydata)]);
eps_y = zeros([1,length(ydata)]);
P_plus = P0;
dx_plus = dx0;

for k = 0:1399

    %dx_vals(:,k+1) = dx; %add perturbation to array

    t = dt*k; %[s] current time

    %Calculate nominal orbit state values at current time
    Xnom = r0*cos(omega0*t);
    Ynom = r0*sin(omega0*t);
    Xnom_dot = -omega0*r0*sin(omega0*t);
    Ynom_dot = omega0*r0*cos(omega0*t);

    xnom = [Xnom;Xnom_dot;Ynom;Ynom_dot]; %nominal state vector at current time

    %true_dx_vals(:,k+1) = x_true(k+1,:)'-xnom; %add state estimate to array

    %store values of LKF stuff
    dX_LKF(:,k+1) = dx_plus;
    X_LKF(:,k+1) = dx_plus + xnom;
    sigma_LKF(:,k+1) = 2*sqrt(diag(P_plus));
    Pk_LKF(:,4*k+1:4*k+4) = P_plus;
    
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
    Gamma_bar = [[0 1 0 0]; [0 0 0 1]]';

    Ft = eye(4)+dt*A_bar; %Ahat(1:4,1:4); %DT state transition matrix at time t
    Omegat = dt*Gamma_bar;
   
    %%%Calculate Linearized Measurements
    Yk = ydata{k+2};
    

    

    if ~isempty(Yk)

        indices = Yk(4,:);
        H = zeros([3*length(indices),4]); %linearized measurement matrix
        Y_vect = zeros([3*length(indices),1]);
        R = zeros([3*length(indices),3*length(indices)]);
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
            R(3*i-2:3*i,3*i-2:3*i) = Rtrue;
            
        end
        
        %%%Time Update Step
        dx_minus = Ft*dx_plus;
        P_minus = Ft*P_plus*Ft' + Omegat*Q*Omegat';
        K = P_minus*H'/(H*P_minus*H' + R);
 
        %dyest(:,k+1) = [H(1:3,1:4)*dx_minus;indices(1)];
    
        %%%Measurement Update Step
        dx_plus = dx_minus + K*(Y_vect-H*dx_minus);
        P_plus = (eye(4) - K*H)*P_minus;
    

        % Calc epsilon_y,k
        eps_y(k+1) = Y_vect'*((H*P_minus*H' + R)\Y_vect);

    else
        %%%Time Update Step
        dx_plus = Ft*dx_plus;
        P_plus = Ft*P_plus*Ft' + Omegat*Q*Omegat';
    end


end


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