%% CONSTANTS
G=6.67384e-11;
M_e=5.97219e24;
R_e=6371000;

rho_SL=1.2;
mean_temp=260;
scale_height=29.26*mean_temp;

%% ROCKET PROPERTIES
CD_roc=0.8;
A_ref=12.6;
m_dry=10000;
m_fuel=120000;
v_exhaust=4000; %or use e.g. 450*9.81

%% ROCKET INITIAL CONDITIONS

%theta=0.0749;
%pos_init=[-1e7,R_e,0];
%vel_init=16000*[cos(theta),sin(theta),0];

pos_init=[0;R_e;0];
vel_init=[1;1;0];

%% ROCKET DESIRED FINAL CONDITIONS
desired_orbenergy=-29600000;

%% SIMULATION PROPERTIES
t_step=0.8;
sim_time=600;

%hillclimbing stuff
max_guesses=500;
missteps_to_reduce=12;
stepsmallifier=1.5;
stpmul=1;
misstepcounter=0;

mdot_schedule=[0,201.593680325808,394.018357219119,30000000;385.801175800646,355.388838869305,89.6945829028151,90];
tvc_schedule=[0,10,86.4027676255286,68.1490413459862,428.296537203446,30000000;0,0.0119371444892641,0.592802205118879,1.14961116538635,1.38462051950640,1.42376917462570]; %currently just t, theta
solution_error=10000;

for hcstep=1:max_guesses
    old_mdot_schedule=mdot_schedule;
    old_tvc_schedule=tvc_schedule;
    old_solution_error=solution_error;
    
    if(hcstep==1)
        kludge=0;
    else
        kludge=1;
    end
    
    mdot_schedule=old_mdot_schedule+kludge*stpmul*[0 10*randn 10*randn 0; 2*randn 2*randn 2*randn 0];
    tvc_schedule=old_tvc_schedule+kludge*stpmul*[0 0 10*randn 10*randn 10*randn 0; 0 0.01*randn 0.01*randn 0.01*randn 0.01*randn 0.01*randn];


    %% RUN SIMULATION
    trajectory=iteratePoweredFlight(pos_init, vel_init, M_e, R_e,  ... % MATLAB, you are utterly vile.
        m_dry, m_fuel, mdot_schedule, tvc_schedule, v_exhaust, CD_roc, A_ref, ...
        rho_SL, scale_height, sim_time, t_step, desired_orbenergy);

    %% CALCULATE ORBITAL ELEMENTS OF FINAL STATE
    orb_elements=orbitalElements(trajectory(2:4,size(trajectory,2)),trajectory(5:7,size(trajectory,2)),M_e);

    %% CALCULATE FITNESS
    solution_error=abs(orb_elements(1)-desired_orbenergy)/3e6+...
        orb_elements(2)*60+...
        (600-trajectory(12,size(trajectory,2)))/600;
    

    
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

    text(100,3000,strcat(num2str(orb_elements(1)),' MJ/kg'));
    text(100,2000,num2str(orb_elements(2)));
    text(100,5000,strcat(num2str(trajectory(12,size(trajectory,2))),' kg'));
    text(100,6000,strcat(num2str(solution_error), ' vs old=', num2str(old_solution_error), ', missteps=', int2str(misstepcounter), ', mult=', num2str(stpmul))); 
    
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
    
    if(solution_error>old_solution_error)
        mdot_schedule=old_mdot_schedule;
        tvc_schedule=old_tvc_schedule;
        solution_error=old_solution_error;
        misstepcounter=misstepcounter+1;
    else
        misstepcounter=0;
    end;
    
    if (misstepcounter>=missteps_to_reduce)
        stpmul=stpmul/stepsmallifier;
        misstepcounter=0;
    end;
end