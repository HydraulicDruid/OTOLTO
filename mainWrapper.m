%% CONSTANTS
G=6.67384e-11;
M_e=5.97219e24;
R_e=6371000;

rho_SL=1.2;
mean_temp=260;
scale_height=29.26*mean_temp;

%% ROCKET PROPERTIES
CD_roc=1.15;
A_ref=1;
m_roc=1000;

%% ROCKET INITIAL CONDITIONS
pos_init=[0,50000+R_e,0];
vel_init=[6000,0,0];

%% SIMULATION PROPERTIES
t_step=1;
sim_time=4;

%% RUN SIMULATION
trajectory=iterateBallisticTrajectory(pos_init, vel_init, M_e, R_e, m_roc,...   % MATLAB, you are utterly vile.
    CD_roc, A_ref, rho_SL, scale_height, sim_time, t_step);