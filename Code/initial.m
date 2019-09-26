clear,clc
tf = 10; %total simulation time
Ts = 0.001;
qi = [0;0;0;0;0;0];% initial configuration
sim UR5Model;
% Animation