clear all;clc;close all;

theta_0 = deg2rad(0); 
theta_f = deg2rad(120);
theta_target = deg2rad(60) ;
vel_target = 10;
max_v = 10;
dt = 0.001;

lfds_plan(theta_0, theta_f, theta_target, vel_target, max_v, dt)
