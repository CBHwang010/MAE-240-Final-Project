%% MAE 240 Final Project - Skeleton Code
% J2 and lunisolar perturbations near the Laplace radius.
% Uses MAE 240 demo functions: koe2rv_local, rv2coe_local, kepEqn_local.

clear; clc; close all;
format compact;

%% Constants

C.mu_E = 398600.435436;      % km^3/s^2
C.R_E  = 6378.1366;          % km
C.J2   = 1.0826359e-3;

C.mu_M = 4902.800066;        % km^3/s^2
C.a_M  = 383397.7725;        % km
C.i_M  = deg2rad(5.145);     % rad

C.mu_S = 1.32712440018e11;   % km^3/s^2
C.a_S  = 149597870.7;        % km
C.eps  = deg2rad(23.439);    % rad

C.n_M = sqrt((C.mu_E + C.mu_M)/C.a_M^3);
C.n_S = sqrt(C.mu_S/C.a_S^3);

%% Adjustable case settings

aRatio = 7.7;        % Change to 6, 7, 7.7, 8, 9, or 10

modelList = {'2B','J2','LS','FULL'};
% 2B   = two-body
% J2   = two-body + J2
% LS   = two-body + Moon + Sun
% FULL = two-body + J2 + Moon + Sun

e0   = 0.01;         % Initial eccentricity
inc0 = deg2rad(30);  % Initial inclination
Om0  = deg2rad(0);   % Initial RAAN
w0   = deg2rad(0);   % Initial argument of perigee
M0   = deg2rad(0);   % Initial mean anomaly

tYears = 5;          % Propagation time
Nt     = 2000;       % Output points

%% Initial state

a0 = aRatio*C.R_E;
coe0 = [a0; e0; inc0; Om0; w0; M0];

% TODO: copy koe2rv_local from mae240_GVE_demo.m
[r0, v0] = koe2rv_local(coe0, C.mu_E);

x0 = [r0; v0];

tFinal = tYears*365.25*86400;
tspan  = linspace(0, tFinal, Nt);

opts = odeset('RelTol',1e-11,'AbsTol',1e-12);

%% Run one semi-major-axis case

Results = struct();

for im = 1:length(modelList)

    modelType = modelList{im};
    fprintf('Running %s model for a/RE = %.2f\n', modelType, aRatio);

    [T, X] = ode45(@(t,x) eom_project(t, x, C, modelType), ...
                   tspan, x0, opts);

    coeHist = zeros(length(T),6);

    for k = 1:length(T)
        r_k = X(k,1:3).';
        v_k = X(k,4:6).';

        % TODO: copy rv2coe_local from mae240_GVE_demo.m
        coeHist(k,:) = rv2coe_local(r_k, v_k, C.mu_E).';
    end

    Results.(modelType).T = T;
    Results.(modelType).X = X;
    Results.(modelType).coeHist = coeHist;

end

%% Plot element changes

figure('Color','w','Name','One Case Comparison');

for im = 1:length(modelList)

    modelType = modelList{im};
    T = Results.(modelType).T;
    coeHist = Results.(modelType).coeHist;

    timeYears = T/(365.25*86400);

    aHist  = coeHist(:,1);
    eHist  = coeHist(:,2);
    iHist  = coeHist(:,3);
    OmHist = unwrap(coeHist(:,4));
    wHist  = unwrap(coeHist(:,5));

    de  = eHist - eHist(1);
    di  = rad2deg(iHist - iHist(1));
    dOm = rad2deg(OmHist - OmHist(1));
    dw  = rad2deg(wHist - wHist(1));

    subplot(2,2,1)
    plot(timeYears, de, 'LineWidth', 1.3); hold on; grid on; box on;
    xlabel('Time [yr]'); ylabel('\Delta e'); title('Eccentricity');

    subplot(2,2,2)
    plot(timeYears, di, 'LineWidth', 1.3); hold on; grid on; box on;
    xlabel('Time [yr]'); ylabel('\Delta i [deg]'); title('Inclination');

    subplot(2,2,3)
    plot(timeYears, dOm, 'LineWidth', 1.3); hold on; grid on; box on;
    xlabel('Time [yr]'); ylabel('\Delta \Omega [deg]'); title('RAAN');

    subplot(2,2,4)
    plot(timeYears, dw, 'LineWidth', 1.3); hold on; grid on; box on;
    xlabel('Time [yr]'); ylabel('\Delta \omega [deg]'); title('Argument of Perigee');

end

for p = 1:4
    subplot(2,2,p)
    legend(modelList,'Location','best');
end

sgtitle(sprintf('Osculating Element Changes for a = %.2f R_E', aRatio));

%% Summary metrics

fprintf('\nSummary for a/RE = %.2f\n', aRatio);
fprintf('----------------------------------------\n');

for im = 1:length(modelList)

    modelType = modelList{im};
    coeHist = Results.(modelType).coeHist;

    aHist  = coeHist(:,1);
    eHist  = coeHist(:,2);
    iHist  = coeHist(:,3);
    OmHist = unwrap(coeHist(:,4));
    wHist  = unwrap(coeHist(:,5));

    da  = aHist - aHist(1);
    de  = eHist - eHist(1);
    di  = rad2deg(iHist - iHist(1));
    dOm = rad2deg(OmHist - OmHist(1));
    dw  = rad2deg(wHist - wHist(1));

    fprintf('\nModel: %s\n', modelType);
    fprintf('  max |Delta a|     = %.4e km\n',  max(abs(da)));
    fprintf('  max |Delta e|     = %.4e\n',     max(abs(de)));
    fprintf('  max |Delta i|     = %.4e deg\n', max(abs(di)));
    fprintf('  max |Delta RAAN|  = %.4e deg\n', max(abs(dOm)));
    fprintf('  max |Delta omega| = %.4e deg\n', max(abs(dw)));

end

%% Local functions

function dxdt = eom_project(t, x, C, modelType)

r = x(1:3);
v = x(4:6);

a = -C.mu_E*r/norm(r)^3;

switch upper(modelType)

    case '2B'

    case 'J2'
        a = a + accel_J2_project(r, C);

    case 'LS'
        rM = moon_position_simple(t, C);
        rS = sun_position_simple(t, C);

        a = a + accel_thirdbody_project(r, rM, C.mu_M);
        a = a + accel_thirdbody_project(r, rS, C.mu_S);

    case 'FULL'
        rM = moon_position_simple(t, C);
        rS = sun_position_simple(t, C);

        a = a + accel_J2_project(r, C);
        a = a + accel_thirdbody_project(r, rM, C.mu_M);
        a = a + accel_thirdbody_project(r, rS, C.mu_S);

    otherwise
        error('Unknown modelType. Use 2B, J2, LS, or FULL.');
end

dxdt = [v; a];

end

function aJ2 = accel_J2_project(r, C)

x = r(1);
y = r(2);
z = r(3);

rmag = norm(r);
r2 = rmag^2;
z2 = z^2;

factor = 3*C.J2*C.mu_E*C.R_E^2/(2*rmag^5);

aJ2 = factor * [x*(5*z2/r2 - 1);
                y*(5*z2/r2 - 1);
                z*(5*z2/r2 - 3)];

end

function a3B = accel_thirdbody_project(rSC, rB, muB)

rho = rB - rSC;
a3B = muB*(rho/norm(rho)^3 - rB/norm(rB)^3);

end

function rM = moon_position_simple(t, C)

theta = C.n_M*t;

rM = C.a_M * [cos(theta);
              sin(theta)*cos(C.i_M);
              sin(theta)*sin(C.i_M)];

end

function rS = sun_position_simple(t, C)

theta = C.n_S*t;

rS = C.a_S * [cos(theta);
              sin(theta)*cos(C.eps);
              sin(theta)*sin(C.eps)];

end

%% TODO: Paste these demo functions below
% From mae240_GVE_demo.m:
%   koe2rv_local
%   rv2coe_local
%   kepEqn_local
%
% Convention:
%   coe = [a; e; i; Omega; omega; M]
%   angles in radians
