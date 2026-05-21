clear; clc; close all;
tau = 8e-9; % Laser pulse FWHM [s]
T0 = 1800; % Initial gas/particle temperature [K]
Dp0 = 16.5e-9; % Initial soot diameter [m]

Lambda_L = 1064e-9; % Laser wavelength [m]
Lambda_d = 450e-9; % Detection wavelength [m]

E_m = 0.3; % Absorption function

% Fluence range [J/cm^2]
F_cm2 = linspace(0.01,0.5,40);
% SI [J/m^2]
F = F_cm2 * 1e4;

t_end = 500e-9;
dt = 0.2e-9;

tvec = 0:dt:t_end;

% Gate timing Bladh
t1 = 2e-9;
gate_width = 20e-9;

[SLII_norm,T,Dp,SLII_tmp,t_out,Tp_max] = Bladh_LII_Model_1( tau,tvec, E_m,T0,Dp0,F,Lambda_L,Lambda_d,gate_width,t1);

figure
plot(F_cm2,SLII_norm,'LineWidth',2)
xlabel('Laser Fluence [J/cm^2]')
ylabel('Normalized Peak LII Signal')
title('Bladh Simplified LII Model')
grid on

figure
hold on
for i = 1:length(F)
    plot(t_out*1e9,T(:,i),'LineWidth',2)
end
xlabel('Time [ns]')
ylabel('Temperature [K]')
title('Particle Temperature Evolution')
legend(string(F_cm2) + " J/cm^2")
grid on

figure
hold on
for i = 1:length(F)
    S = SLII_tmp(:,i);
    S = S ./ max(S);
    plot(t_out*1e9,S,'LineWidth',2)
end
xlabel('Time [ns]')
ylabel('Normalized LII Signal')
title('Temporal LII Signal')
legend(string(F_cm2) + " J/cm^2")
xlim([0 500])
grid on

figure
hold on
for i = 1:length(F)
    plot(t_out*1e9,Dp(:,i)*1e9,'LineWidth',2)
end
xlabel('Time [ns]')
ylabel('Particle Diameter [nm]')
title('Soot Diameter Evolution')
legend(string(F_cm2) + " J/cm^2")
grid on

figure
plot(F_cm2,Tp_max,'LineWidth',2)
xlabel('Laser Fluence [J/cm^2]')
ylabel('Peak Particle Temperature [K]')
title('Peak Temperature vs Fluence')
grid on

function [SLII_norm,T,Dp,SLII_tmp,t_out,Tp_max] = ...
    Bladh_LII_Model_1( ...
    tau,tvec,E_m,T0,Dp0,F,...
    Lambda_L,Lambda_d,dt,t1)
    dtau = 1.5*tau;
    h = 6.626e-34;
    c = 2.997e8;
    kB = 1.38e-23;

    for i = 1:numel(F)
        I_peak = F(i)*sqrt((4*log(2))/(pi*tau^2));

        q_t = @(t) I_peak .* exp((-4.*log(2).*((t-dtau)./tau).^2));
        X0 = [T0; Dp0];
        opts = odeset( 'RelTol',1e-8,'AbsTol',1e-9, 'MaxStep',0.1e-9);

        [t_out,X] = ode45( @(t,X) LII_Model(t,X,q_t,Lambda_L,E_m,T0),tvec,X0,opts);

        T(:,i) = X(:,1);
        Dp(:,i) = X(:,2);

        t2 = t1 + dt;

        [~,idx1] = min(abs(t_out-t1));
        [~,idx2] = min(abs(t_out-t2));
        Tp_max(i) = max(T(:,i));

        SLII_tmp(:,i) = (8.*pi^3.*Dp(:,i).^3.*h.*c.^2.*E_m) ./ (Lambda_d.^6 .*  (exp(h.*c ./ (Lambda_d.*kB.*T(:,i))) - 1));
        SLII_int(i) =  trapz( t_out(idx1:idx2),SLII_tmp(idx1:idx2,i))  /(t2-t1);
    end
    SLII_norm = SLII_int ./ max(SLII_int);

end

function dXdt = ...
    LII_Model(t,X,q_t,lambda,E_m,Tg)

    T = X(1);
    Dp = X(2);
    P0 = 101325;
    R = 8.314;
    h = 6.626e-34;
    c = 2.997e8;
    kB = 1.38e-23;
    rho = 2200;
    Cs = 2100;
    Vp = (pi/6)*Dp^3;
    Ap = pi*Dp^2;
    mp = rho*Vp;
    q = q_t(t);

    Qabs = (pi^2*Dp^3*E_m/lambda)* q;

    ka = 0.026;
    lambda_MFP = 0.5665e-6;
    G = 22.064;
    Qcond = 2*ka*Ap*(T-Tg) /(Dp + G*lambda_MFP);

    Qrad = ((199*pi^3*Dp^3*E_m) /(h*(h*c)^3)) * ((kB*T)^5 - (kB*Tg)^5);

    Pv = ...
        exp( ...
        -122.96 ...
        + 9.0558e-2*T ...
        - 2.7637e-5*T^2 ...
        + 4.1754e-9*T^3 ...
        - 2.4875e-13*T^4) ...
        * 101325;

    Mv = ...
        (17.179 ...
        + 6.8654e-4*T ...
        + 2.9962e-6*T^2 ...
        - 8.5954e-10*T^3 ...
        + 1.0486e-13*T^4) ...
        /1000;

    dHv = ...
        2.05398e5 ...
        + 7.3660e2*T ...
        - 0.40713*T^2 ...
        + 1.1992e-4*T^3 ...
        - 1.7946e-8*T^4 ...
        + 1.0717e-12*T^5;

    alpha_m = 0.8;

    dMdt =  pi*Dp^2 * alpha_m  * Pv  * sqrt(Mv/(2*pi*R*T));
    Qsub = dHv*dMdt/Mv;
    dTdt = (Qabs - Qcond - Qrad - Qsub)/(mp*Cs);

    dDpdt = -2*dMdt /(rho*pi*Dp^2);
    dXdt = [dTdt; dDpdt];

end