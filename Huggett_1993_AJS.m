% Huggett 1993 Replication
% Macro- ECON 516 Midterm
% Aditi Singh, Jan Rosa, Sudipta Ghosh 


tic
clear all

%% 1. Parameters

% endowment
eh=1; 
el=0.1;
pihh=0.925;
pihl=0.5;
trans_mat=[pihh pihl; 1-pihh 1-pihl]; % this is the transition matrix [phh phl;plh pll]

beta=0.99322;
sigma=1.5; % risk aversion parameter
sigma1=3;
% asset grid
amin_c=[-2 -4 -6 -8];
amax=8;
grid_len=500; % # of grid points 

% tolarence levels
c_tol=1e-7;
q_tol=2.5e-4;

%% Try function for different sigmas

[r_15, q_15] = HuggettSolve(amin_c,amax,grid_len,eh,el,beta, sigma,c_tol,q_tol,trans_mat);
[r_3, q_3] = HuggettSolve(amin_c,amax,grid_len,eh,el,beta, sigma1,c_tol,q_tol,trans_mat);

 
