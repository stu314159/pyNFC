% File Processing Script for Subspace data
clear;
clc;
close all;

% load the data:
load('ss_data.mat');

% the environment will now contain ss_rho, ss_ux, ss_uy, and ss_uz
% get some parameters
[Num_ts, Nz, Ny] = size(ss_rho);

ux = nan(Nz,Ny,Num_ts);
uy = nan(Nz,Ny,Num_ts);
uz = nan(Nz,Ny,Num_ts);
pressure = nan(Nz,Ny,Num_ts);

for ts = 1:Num_ts
  % attempt to reorder the data
  ux(:,:,ts) = ss_ux(ts,:,:);
  uy(:,:,ts) = ss_uy(ts,:,:);  
  uz(:,:,ts) = ss_uz(ts,:,:);
  pressure(:,:,ts) = ss_rho(ts,:,:);
end

% make plots of uz (for example)
for ts = 1:Num_ts
  imagesc(uz(:,:,ts))
  axis equal
  pause(1)
  close 'all'
end
