% iterateBallisticTrajectory(position, velocity, priMass, priRad, 
%   m_obj, C_D, A_ref, rho_SL, scaleheight,
%   simTime, timestep)

% Note iterationTime = 0 gives one-period case and iterationTime = -1
% gives until-ground-hit case.

% returns a matrix whose columns are timesteps and whose eleven rows are:
% timestep
% position[3]
% velocity[3]
% acceleration[3]
% Q
% 
% things to be added later:
% something about atmospheric heating

function trajectory = iterateBallisticTrajectory(pos_init, vel_init, priMass, priRad, m_obj, C_D, A_ref, rho_SL, scaleheight, simTime, timestep)

R_min=1000;

t=0;
stepsmade=1;
pos=pos_init;
vel=vel_init;

trajectory=zeros(11,1);
trajectory(2:4,1)=pos;
trajectory(5:7,1)=vel;

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
    
    %Use instantaneous acceleration to decide on timestep
    dt=(timestep*9.81)/norm(pwrd_acc(pos,vel,m_obj,priMass,C_D,A_ref,air_density,0));
    trajectory(1,stepsmade)=dt;   
    trajectory(8:10,stepsmade)=pwrd_acc(pos,vel,m_obj,priMass,C_D,A_ref,air_density,0);

    %Pick up position and velocity for this timestep
    pos=trajectory(2:4,stepsmade);
    vel=trajectory(5:7,stepsmade);
    
    %Then work out Q, because things are all done in a stupid order anyway
    trajectory(11,stepsmade)=1/2*norm(vel)*norm(vel)*air_density;
    
    %RK4 constants
    kv1=dt*pwrd_acc(pos,vel,m_obj,priMass,C_D,A_ref,air_density,0);
    kx1=dt*vel;
    
    kv2=dt*pwrd_acc(pos+kx1/2,vel+kv1/2,m_obj,priMass,C_D,A_ref,air_density,0);
    kx2=dt*(vel+kv1/2);
    
    kv3=dt*pwrd_acc(pos+kx2/2,vel+kv2/2,m_obj,priMass,C_D,A_ref,air_density,0);
    kx3=dt*(vel+kv2/2);

    kv4=dt*pwrd_acc(pos+kx3,vel+kv3,m_obj,priMass,C_D,A_ref,air_density,0);
    kx4=dt*(vel+kv3);
    
    %Log next values of position and velocity.
    trajectory(2:4,stepsmade+1)=pos+(kx1+2*kx2+2*kx3+kx4)/6;
    trajectory(5:7,stepsmade+1)=vel+(kv1+2*kv2+2*kv3+kv4)/6;

    t=t+dt;
    stepsmade=stepsmade+1;    
    
end;
