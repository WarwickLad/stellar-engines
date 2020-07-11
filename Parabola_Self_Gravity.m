%% clear contents

clear all
close all
clc

%% constants %%

pkg load symbolic
format long

au = 1.495978707E11;
G = 6.67408E-11;
Ag_density = 10.5E3;
Fe_density = 7.8E3;
Al_density = 2.7E3;
thickness = 0.2E-6;

F = 0.4.*au; % focus of parabola
area_total = 4*pi*F;
constant = -1.5*thickness*Al_density*G;

%% Specify fine-grainness %%

for N = [11]; % course-grainness of simulation surface

  %% Plot Paraboloid %%

  theta = linspace(pi/2,3*pi/2,N); % azimuthal angle
  phi = linspace(0,2*pi,N); % polar angle
  r = (2.*F)./(1-cos(theta)); % paraboloid surface
  [THETA,PHI] = meshgrid(theta,phi);

  R = (2.*F)./(1-cos(THETA));
  X = R.*sin(THETA).*(cos(PHI));
  Y = R.*sin(THETA).*(sin(PHI));
  Z = R.*cos(THETA);

  %% Calculate area units %%
  
  d_theta = pi/(N-1);
  d_phi = 2*pi/(N-1);
  Area = R.^2.*d_theta.*d_phi;
  
  %% Calculate Gravity %%

  self_gravity_x = zeros(N,N);
  self_gravity_y = zeros(N,N);
  self_gravity_z = zeros(N,N);

  for i = [1:N];
    for j = [1:N];
      
      function_Z_12 = zeros(1,N);
      
      X1 = X - X(i,j);
      Y1 = Y - Y(i,j);
      Z1 = Z - Z(i,j);
      
      k = -1.5;
      
      RHO = ((X1.^2 + Y1.^2 + Z1.^2).^(k));
      RHO(RHO > 1E10) = 0;
      
      function_X_1 = Area(i,j).*Area.*constant.*RHO.*X1;
      function_Y_1 = Area(i,j).*Area.*constant.*RHO.*Y1;
      function_Z_1 = Area(i,j).*Area.*constant.*RHO.*Z1;
      
      self_gravity_x(i,j) = sum(function_X_1(:));
      self_gravity_y(i,j) = sum(function_Y_1(:));
      self_gravity_z(i,j) = sum(function_Z_1(:));
      
    endfor
    
  endfor
  
  for i = [1:N];
    self_gravity_z(:,i) = mean(self_gravity_z(:,i));
  endfor
  
  self_gravity = (self_gravity_x.^2 + self_gravity_y.^2 + self_gravity_z.^2).^(0.5);
  self_acceleration = self_gravity./(thickness.*Al_density./Area);
  self_acceleration_x_abs = abs(self_gravity_x)./(thickness.*Al_density./Area);
  self_acceleration_y_abs = abs(self_gravity_y)./(thickness.*Al_density./Area);
  self_acceleration_z_abs = abs(self_gravity_z)./(thickness.*Al_density./Area);
  
  max(self_acceleration_x_abs(:))
  max(self_acceleration_y_abs(:))
  max(self_acceleration_z_abs(:))

  figure; plot(Z,self_acceleration_z_abs,'bx')
  
endfor