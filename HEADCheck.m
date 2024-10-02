function [Power] = HEADCheck(laJ211, aaJ211, fs, sex)
% -------------------------------------------------------------------------
% Version: 1.0 (Released: Oct 2024)
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
% Research contact: Dr Gregory Tierney (g.tierney@ulster.ac.uk)
% -------------------------------------------------------------------------
% Input variable:
% laJ211: Componential linear acceleration (X,Y,Z) in J211 coordinate system
% aaJ211: Componential angular acceleration (X,Y,Z) in J211 coordinate system
% fs: Sampling frequency
% sex: Enter 'm' for male or 'f' for female
% -------------------------------------------------------------------------
% Output variable: 
% Power: Rate of change of kinetic energy that the head undergoes during a HAE
% -------------------------------------------------------------------------

% Time Parameters
T = 1/fs; % Period (s)
L = length(laJ211); % Data length
Time = (0:L-1)*T'; % Time vector (s)
dT = Time(2)-Time(1); % Time step (s)

% Mass
if sex=='m'; mass = 4.1; elseif sex=='f'; mass = 3.2; end 

% Moments of Inertia  
Ixx= (74.8 * mass -125.5)/(100^2); Iyy= (71.4 * mass -90.2)/(100^2); Izz= (45.6 * mass -26.5)/(100^2);

% Power 
Power = (mass.*laJ211(:,1).*(cumtrapz(laJ211(:,1))*dT)) + (mass.*laJ211(:,2).*(cumtrapz(laJ211(:,2))*dT)) + (mass.*laJ211(:,3).*(cumtrapz(laJ211(:,3))*dT)) + (Ixx.*aaJ211(:,1).*(cumtrapz(aaJ211(:,1))*dT)) + (Iyy.*aaJ211(:,2).*(cumtrapz(aaJ211(:,2))*dT)) + (Izz.*aaJ211(:,3).*(cumtrapz(aaJ211(:,3))*dT));