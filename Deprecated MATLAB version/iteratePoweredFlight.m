% iterateBallisticTrajectory(position, velocity, priMass, priRad, 
%   m_dry, m_fuel, thrustschedule, tvcschedule, v_exh, C_D, A_ref, 
%   rho_SL, scaleheight,
%   simTime, timestep, stopatmeco)

% Note iterationTime = 0 gives one-period case and iterationTime = -1
% gives until-ground-hit case.

% returns a matrix whose columns are timesteps and whose twelve rows are:
% timestep
% position[3]
% velocity[3]
% acceleration[3]
% Q
% remaining propellant
% 
% things to be added later:
% something about atmospheric heating

function trajectory = iteratePoweredFlight(pos_init, vel_init, priMass, priRad, m_dry, m_fuel, thrustschedule, tvcschedule, v_exh, C_D, A_ref, rho_SL, scaleheight, simTime, timestep, desired_orbenergy, stopatmeco)

R_min=priRad-100;

t=0;
stepsmade=1;
pos=pos_init;
vel=vel_init;

meco=0;

trajectory=zeros(12,1);
trajectory(2:4,1)=pos;
trajectory(5:7,1)=vel;
trajectory(12,1)=m_fuel;

m_obj=m_dry+m_fuel;

while(t<simTime)
    r=norm(pos);
    if(r<R_min)
        break;
    end;
    
    if(r>priRad)
        air_density=rho_SL*exp(-(r-priRad)/scaleheight);
    else
        air_density=0;
    end;
    
    orb_elements=orbitalElements(trajectory(2:4,stepsmade),trajectory(5:7,stepsmade),priMass);

    if(orb_elements(1)>desired_orbenergy)
        if(stopatmeco)
            break;
        end
        meco=1;    
    end;
    
    thrusttheta=interp1(tvcschedule(1,:),tvcschedule(2,:),t);
    thr_uvec=[sin(thrusttheta);cos(thrusttheta);0];
    
    %Use instantaneous acceleration to decide on timestep
    thrust=emdot(t,thrustschedule,m_obj,m_dry,meco)*v_exh*thr_uvec;
    
    dt=(timestep*9.81)/norm(pwrd_acc(pos,vel,m_obj,priMass,C_D,A_ref,air_density,thrust));
    trajectory(1,stepsmade)=dt;   
    trajectory(8:10,stepsmade)=pwrd_acc(pos,vel,m_obj,priMass,C_D,A_ref,air_density,thrust);

    %Pick up position and velocity for this timestep
    pos=trajectory(2:4,stepsmade);
    vel=trajectory(5:7,stepsmade);
    
    %Then work out Q, because things are all done in a stupid order anyway
    
    trajectory(11,stepsmade)=1/2*norm(vel)*norm(vel)*air_density*norm(vel);
    
    %RK4 constants#
    km1=dt*emdot(t,thrustschedule,m_obj,m_dry,meco);
    thrustacc1=v_exh*km1/dt;
    thrusttheta=interp1(tvcschedule(1,:),tvcschedule(2,:),t);
    thr_uvec=[sin(thrusttheta);cos(thrusttheta);0];
    kv1=dt*pwrd_acc(pos,vel,m_obj,priMass,C_D,A_ref,air_density,thrustacc1*thr_uvec);
    kx1=dt*vel;
    
    km2=dt*emdot(t+dt/2,thrustschedule,m_obj-(km1/2),m_dry,meco);
    thrustacc2=v_exh*km2/dt;    
    thrusttheta=interp1(tvcschedule(1,:),tvcschedule(2,:),t+dt/2);
    thr_uvec=[sin(thrusttheta);cos(thrusttheta);0];
    kv2=dt*pwrd_acc(pos+kx1/2,vel+kv1/2,m_obj,priMass,C_D,A_ref,air_density,thrustacc2*thr_uvec);
    kx2=dt*(vel+kv1/2);
    
    km3=dt*emdot(t+dt/2,thrustschedule,m_obj-(km2/2),m_dry,meco);
    thrustacc3=v_exh*km3/dt;
    thrusttheta=interp1(tvcschedule(1,:),tvcschedule(2,:),t+dt/2);
    thr_uvec=[sin(thrusttheta);cos(thrusttheta);0];
    kv3=dt*pwrd_acc(pos+kx2/2,vel+kv2/2,m_obj,priMass,C_D,A_ref,air_density,thrustacc3*thr_uvec);
    kx3=dt*(vel+kv2/2);

    km4=dt*emdot(t+dt,thrustschedule,m_obj-km3,m_dry,meco);
    thrustacc4=v_exh*km3/dt;    
    thrusttheta=interp1(tvcschedule(1,:),tvcschedule(2,:),t+dt);
    thr_uvec=[sin(thrusttheta);cos(thrusttheta);0];    
    kv4=dt*pwrd_acc(pos+kx3,vel+kv3,m_obj,priMass,C_D,A_ref,air_density,thrustacc4*thr_uvec);
    kx4=dt*(vel+kv3);
    
    %Log next values of position and velocity.
    trajectory(2:4,stepsmade+1)=pos+(kx1+2*kx2+2*kx3+kx4)/6;
    trajectory(5:7,stepsmade+1)=vel+(kv1+2*kv2+2*kv3+kv4)/6;
    m_obj=m_obj-(km1+2*km2+2*km3+km4)/6;
    trajectory(12,stepsmade+1)=m_obj-m_dry;

    t=t+dt;
    stepsmade=stepsmade+1;    
    
end;
