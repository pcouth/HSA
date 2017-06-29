%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% A script to test HSA generation for typical 1D1P, 1DMP, MD1P, MDMP cases
% using various options
%--------------------------------------------------------------------------
close all; clear all; clc

%% One-Dimension, One-Parameter (1D1P)

% Saddle Node Bifurcation
f1 = @(x,u,p) p+x^2;  
x01 = -1; 
p01 = -1;
c1 = cosy(f1,1,0,1);
c1.p_min = -1;              
c1.p_max = 1;
[b1,bp1] = c1.trace_solution_branches(x01,[],p01);

node1fp = hsa1D('SN_fp',f1,b1,bp1,-0.5,'dot2tex');           
node1gen = hsa1D('SN_gen',f1,b1,bp1,[],'dot2tex');           
node1genR = hsa1D('SN_gen',f1,b1,bp1,[],{'x' '\gamma ';'p' '\alpha '});

% Pitchfork Bifurcation
f2 = @(x,u,p) p*x-x^3; 
x02 = 0; 
p02 = -1;
c2 = cosy(f2,1,0,1);
c2.p_min = -1;              
c2.p_max = 1;
[b2,bp2] = c2.trace_solution_branches(x02,[],p02);
node2fp = hsa1D('PF_fp',f2,b2,bp2,0.5);             
node2gen = hsa1D('PF_gen',f2,b2,bp2,[]);            

% Transcritical Bifurcation
f3 = @(x,u,p) p*x-x^2; 
x03 = 0; 
p03 = -1;
c3 = cosy(f3,1,0,1);
c3.p_min = -1;              
c3.p_max = 1;
[b3,bp3] = c3.trace_solution_branches(x03,[],p03);
node3fp = hsa1D('TC_fp',f3,b3,bp3,-0.5);            
node3gen = hsa1D('TC_gen',f3,b3,bp3,[]);            

%% One-Dimension, Multi-Parameter (1DMP)

f4 = @(x,u,p) p(1)+p(2)*x+p(3)*x^2-p(4)*x^3;
PD4 = [ NaN 0 0 -1 ; 0 NaN 0 1 ; 0 NaN 1 1 ];
x04 = [-1;-1;-1]; 
p04 = [-1;-1;-1];
plim4 = [-1 1;-1 1;-1 1];
[F4,B4,BP4] = dataMP(f4,PD4,x04,p04,plim4);    
opt4.PD = PD4;
opt4.dot2tex = 'on';
opt4.texlab = {'x' '\gamma';'p_1' '\alpha';'p_2' 'v';'p_3' 'h';'p_4' '\rho'};
% opt4.ztol = 0.0001;
node4fp = hsa1D('ODMP_fp',F4,B4,BP4,0.5,PD4);   	
node4gen = hsa1D('ODMP_gen',F4,B4,BP4,[],opt4);  	

%% Multi-Dimension, One-Parameter (MD1P)

f5 = @(x,u,p) [ -p*x(2)-x(1)^3 ; -x(2)^2-x(1) ];
x05 = [-1;-1]; 
p05 = -1;
c5 = cosy(f5,2,0,1);
c5.p_min = -10; 
c5.p_max = 10;
c5.x_max = [10;10]; 
c5.x_min = [-10;-10];
c5.cm.max_steps = 15000;
[b5,bp5] = c5.trace_solution_branches(x05,[],p05);
node5fp = hsaMD('MD1P_fp',f5,b5,bp5,-0.5);       	
node5gen = hsaMD('MD1P_gen',f5,b5,bp5,[],'dot2tex');       	
node5genR = hsaMD('MD1P_gen',f5,b5,bp5,[],{'x_1' '\gamma ';'x_2' 'v';'p' '\alpha '});

%% Multi-Dimension, Multi-Parameter (MDMP)

f6 = @(x,u,p) [ -p(1)*x(2)+p(2)*x(1)-x(1)^3 ; -p(3)*(x(2)^2+x(1))+p(2)*x(1)-p(4)*x(1)^2 ];
PD6 = [ NaN 0 1 0 ; 1 NaN 0 1 ];
x06 = [0 0;0 0];
p06 = [-1;-1];
plim6 = [-10 10;-10 10];
[F6,B6,BP6] = dataMP(f6,PD6,x06,p06,plim6);  
opt6.PD = PD6;
opt6.dot2tex = 'on';
node6fp = hsaMD('MDMP_fp',F6,B6,BP6,-0.5,opt6);      
opt6.texlab = {'x_1' '\gamma ';'x_2' 'v'};
node6gen = hsaMD('MDMP_gen',F6,B6,BP6,[],opt6);    	

xlimit = [-3 3;-3 3];
tspan = [0 20];
segs = 10; 
basin6 = basinMD(node6fp,'Simulation',xlimit,segs,tspan);


