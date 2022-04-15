% Gives thermodynamic contants for given air temperature (temp) in Celcius
% Based on “Acoustical wave propagation in cylindrical ducts: Transmission 
% line parameter approximations for isothermal and nonisothermal boundary 
% conditions”, D. Keefe, J. Acoust. Soc. Am. 75, 58–62 (1984).
%@param temp Temperature in Celsius
%@return vector of thermodynamic constants [air density, speed of sound,
% shear viscosity, Prandtl number, ratio of specific heats]
function [rho, c, eta, nu, gamma] = ThermConstants(temp)
deltaT=temp-26.85;%temperature difference
rho=1.1769*(1-0.00335*deltaT);%air density
c=347.23*(1+0.00166*deltaT);%speed of sound
eta=1.846e-5*(1+0.0025*deltaT);%shear viscosity
nu=0.841*(1-0.0002*deltaT);%prandtl
gamma=1.4017*(1-0.00002*deltaT);%ratio of specific heats
