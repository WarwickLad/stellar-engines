%% clear contents

clear all
close all
clc

%% constants %%

pkg load symbolic
pkg load parallel
format long

au = 1.495978707E11;
G = 6.67408E-11;
Ag_density = 10.5E3;
Fe_density = 7.8E3;
Al_density = 2.7E3;
thickness = 0.2E-6;

self_acceleration_x_abs_max = zeros(1,5);
self_acceleration_y_abs_max = zeros(1,5);
self_acceleration_z_abs_max = zeros(1,5);
##F = [10,0.2.*au,0.4.*au,0.6.*au,0.8.*au];
F = [0.2.*au,0.4.*au];

for ij = [1:length(F)]; % focus of parabola
  
  area_total = 4*pi*F(ij);
  constant = 1.5;

  %% Specify fine-grainness %%
  
  N = 127; % course-grainness of simulation surface

  %% Paraboloid Areas %%

  N = N + 1;
  
  theta = linspace(pi/2,3*pi/2,N); % azimuthal angle
  phi = linspace(0,2*pi,N); % polar angle
  r = (2.*F(ij))./(1-cos(theta)); % paraboloid surface
  [THETA,PHI] = meshgrid(theta,phi);

  R = (2.*F(ij))./(1-cos(THETA));
  X = R.*sin(THETA).*(cos(PHI));
  Z = R.*cos(THETA);

  %% Calculate area units %%
  
  x = X((N/2+1),:);
  z = Z((N/2+1),:);
  
  Area = ones(N-1,N-1);
  
  for i = [1:N/2 - 1];
    Area(i,:) = pi.*(x(i+1)-x(i)).*sqrt(2*(x(i+1)-x(i)).^2+(z(i+1)-z(i)).^2);
    Area(N-i,:) = Area(i,:);
  endfor
  Area(N/2,:) = pi.*(x(N/2+1)-x(N/2)).*sqrt(2*(x(N/2+1)-x(N/2)).^2+(z(N/2+1)-z(N/2)).^2);
  
  %% Paraboloid centre coordinates %%

  N = N - 1;
  
  theta = linspace(pi/2,3*pi/2,N); % azimuthal angle
  phi = linspace(0,2*pi,N); % polar angle
  r = (2.*F(ij))./(1-cos(theta)); % paraboloid surface
  [THETA,PHI] = meshgrid(theta,phi);

  R = (2.*F(ij))./(1-cos(THETA));
  X = R.*sin(THETA).*(cos(PHI));
  Y = R.*sin(THETA).*(sin(PHI));
  Z = R.*cos(THETA);
  
  %% Calculate Gravity %%

  self_gravity_x = zeros(N,N);
  self_gravity_y = zeros(N,N);
  self_gravity_z = zeros(N,N);
  
  for i = [1:N];
    parfor j = [1:N];
      
      function_Z_12 = zeros(1,N);
      
      X1 = X - X(i,j);
      Y1 = Y - Y(i,j);
      Z1 = Z - Z(i,j);
      
      k = -1.5;
      
      RHO = ((X1.^2 + Y1.^2 + Z1.^2).^(k));
      RHO(RHO > 1.5E10) = 0;
      
      function_X_1 = -1.*constant.*Area(i,j).*thickness.*Al_density.*Area.*thickness.*Al_density.*G.*RHO.*X1; % gravity contribution from each paraboloid area in x-direction
      function_Y_1 = -1.*constant.*Area(i,j).*thickness.*Al_density.*Area.*thickness.*Al_density.*G.*RHO.*Y1; % gravity contribution from each paraboloid area in y-direction
      function_Z_1 = -1.*constant.*Area(i,j).*thickness.*Al_density.*Area.*thickness.*Al_density.*G.*RHO.*Z1; % gravity contribution from each paraboloid area in z-direction
      
      function_X_1(i,j) = 0;
      function_Y_1(i,j) = 0;
      function_Z_1(i,j) = 0;
      
      self_gravity_x(i,j) = sum(function_X_1(:)); % total gravity at point x,y,z is sum of x components
      self_gravity_y(i,j) = sum(function_Y_1(:)); % total gravity at point x,y,z is sum of y components
      self_gravity_z(i,j) = sum(function_Z_1(:)); % total gravity at point x,y,z is sum of z components
      
    endparfor
    
  endfor
  
  parfor i = [1:N];
    self_gravity_z(:,i) = mean(self_gravity_z(:,i)); % average small float differences in z-coordinate calculations
  endparfor
  
  self_gravity = (self_gravity_x.^2 + self_gravity_y.^2 + self_gravity_z.^2).^(0.5);
  self_acceleration = self_gravity./(Area.*thickness.*Al_density);
  
  self_acceleration_x_abs = abs(self_gravity_x)./(Area.*thickness.*Al_density);
  self_acceleration_y_abs = abs(self_gravity_y)./(Area.*thickness.*Al_density);
  self_acceleration_z_abs = abs(self_gravity_z)./(Area.*thickness.*Al_density);
  
  self_acceleration_x_abs_max(ij) = max(self_acceleration_x_abs(:));
  self_acceleration_y_abs_max(ij) = max(self_acceleration_y_abs(:));
  self_acceleration_z_abs_max(ij) = max(self_acceleration_z_abs(:));
  
endfor

##figure; surf(X,Y,Z,self_gravity_x); colormap('jet'); shading interp; colorbar; axis image
##figure; surf(Y,X,Z,self_gravity_x); colormap('jet'); shading interp; colorbar; axis image
figure; surf(X,Y,Z,self_gravity_z); colormap('jet'); shading interp; colorbar; axis image

self_acceleration_x_abs_max
self_acceleration_y_abs_max
self_acceleration_z_abs_max

##save self_acceleration_x_max_different_F.mat self_acceleration_x_abs_max
##save self_acceleration_y_max_different_F.mat self_acceleration_y_abs_max
##save self_acceleration_z_max_different_F.mat self_acceleration_z_abs_max