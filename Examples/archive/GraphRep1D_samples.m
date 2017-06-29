%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
% Runs graph1D.m to produce graph representations of 3 prototypical cases

% Note: requires GraphViz, and the addition of the folder containing 
% GraphViz "dot.exe" to environmental variables

% Use 'p = []' for general cases, or specify for fixed-p cases

clear all; close all; clc

%% Saddle Node Bifurcation: x_dot = p+x^2
load('SaddleNode.mat')
label = 'SaddleNode';
p = [];         
graph1D(f,b,bp,p,label)     

%% Pitchfork Bifurcation: x_dot = p*x-x^3
load('Pitchfork.mat')
label = 'Pitchfork';
p = [];  
graph1D(f,b,bp,p,label)      

%% Transcritical Bifurcation: x_dot = p*x-x^2
load('Transcritical.mat')
label = 'Transcritical';
p = [];  
graph1D(f,b,bp,p,label)      
