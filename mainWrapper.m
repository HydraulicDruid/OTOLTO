%% CONSTANTS
G=6.67384e-11;
M_e=5.97219e24;
R_e=6371000;

rho_SL=1.2;
mean_temp=260;
scale_height=29.26*mean_temp;

%% ROCKET PROPERTIES
CD_roc=1.15;
A_ref=0.2;
m_dry=800;
m_fuel=9000;
v_exhaust=4000; %or use e.g. 450*9.81
mdot_schedule=[0 200 400 3e7; 30 30 10 10];
tvc_schedule=[0 10 100 350 3e7; 0 0 pi/6 pi/2 pi/2]; %currently just t, theta

%% ROCKET INITIAL CONDITIONS

%theta=0.0749;
%pos_init=[-1e7,R_e,0];
%vel_init=16000*[cos(theta),sin(theta),0];

pos_init=[0;R_e;0];
vel_init=[1;1;0];

%% ROCKET DESIRED FINAL CONDITIONS
desired_orbenergy=-29600000;

%% SIMULATION PROPERTIES
t_step=1;
sim_time=600;

%% RUN SIMULATION
trajectory=iteratePoweredFlight(pos_init, vel_init, M_e, R_e,  ... % MATLAB, you are utterly vile.
    m_dry, m_fuel, mdot_schedule, tvc_schedule, v_exhaust, CD_roc, A_ref, ...
    rho_SL, scale_height, sim_time, t_step, desired_orbenergy);

%% CALCULATE ORBITAL ELEMENTS OF FINAL STATE
orb_elements=orbitalElements(trajectory(2:4,size(trajectory,2)),trajectory(5:7,size(trajectory,2)),M_e);

%% PLOT SOME PLOTS OR SOMETHING

scrsize=get(0,'ScreenSize');

if (ishandle(trajfig)==false)
    trajfig=figure('OuterPosition',[0 scrsize(4)/2 scrsize(3)/2 scrsize(4)/2]);
    velfig=figure('OuterPosition',[scrsize(3)/2 0 scrsize(3)/2 scrsize(4)/2]);
    accfig=figure('OuterPosition',[0 0 scrsize(3)/2 scrsize(4)/2]);
    qfig=figure('OuterPosition',[scrsize(3)/2 scrsize(4)/2 scrsize(3)/2 scrsize(4)/2]);
end;

figure(trajfig);
title('Trajectory');
clf;
theta=0:pi/50:2*pi;
x=R_e*cos(theta);
y=R_e*sin(theta);

xmax=max(trajectory(2,:));
xmin=min(trajectory(2,:));
ymax=max(trajectory(3,:));
ymin=min(trajectory(3,:));

xcentre=(xmin+xmax)/2;
ycentre=(ymin+ymax)/2;

xmax=xmax+(xmax-xcentre)*0.1;
xmin=xmin-(xcentre-xmin)*0.1;

ymax=ymax+(ymax-ycentre)*0.1;
ymin=ymin-(ycentre-ymin)*0.1;

axis equal;
axis([xmin xmax ymin ymax]);

hold on;
plot(trajectory(2,:),trajectory(3,:),'r');
plot(x,y);
hold off;

%---

figure(qfig);
title('Q');
trajcumsum=cumsum(trajectory,2);
plot(trajcumsum(1,:),trajectory(11,:),'g');

%---

figure(velfig);
title('Velocity');
vels=trajectory(5:7,:);
plot(trajcumsum(1,:),sqrt(sum(vels.^2,1)),'b');

%---
figure(accfig);
title('Acceleration');
accs=trajectory(8:10,:);
plot(trajcumsum(1,:),sqrt(sum(accs.^2,1)),'k');