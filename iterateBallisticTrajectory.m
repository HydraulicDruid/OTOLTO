% iterateBallisticTrajectory(position, velocity, priMass, 
%   objMass, objC_D, objFrontalArea, rho_SL, scaleheight,
%   iterationTime)

% Note iterationTime = 0 gives one-period case and iterationTime = -1
% gives until-ground-hit case

% returns a matrix whose columns are timesteps and whose twelve rows are:
% timestep size
% position[3]
% velocity[3]
% acceleration[3]
% Q
% compression heatflux

function trajectory = iterateBallisticTrajectory(pos, vel, priMass)