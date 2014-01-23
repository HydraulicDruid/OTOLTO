from math import *
import numpy as np
import scipy.interpolate as spip
from scipy.integrate import ode,odeint
import copy
import warnings

G=6.67384e-11
crash_alt=-5000 #no collision detection, so just end the simulation if >5km below ground

def unit_vector(vector):
    """ Returns the unit vector of the vector.
    Code stolen from StackOverflow: 
    http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python  """
    if np.linalg.norm(vector)==0.0:
        return np.zeros(np.shape(vector))
    else:
        return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
            
    Code stolen from StackOverflow: 
    http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

class Atmosphere(object):
    """Currently "type" does nothing; eventually it'll allow switching between
    scale-height and International Standard Atmosphere"""
    def __init__(self,type, rho_SL, scale_height_constant, mean_temp):
        self.type=0
        self.rho_SL=rho_SL
        self.scaleHeight=scale_height_constant*mean_temp
        
    def get_rho(self,mode,altitude):
        """returns either -1 if below ground level or density at altitude if 
        above. Again, mode does nothing, but will eventually allow choice 
        between scale-height model and ISA."""
        if altitude<0:
            return 0
        else:
            return self.rho_SL*exp(-altitude/self.scaleHeight)

class Planet(object):
    """A spherical planet with a mass, radius, atmosphere, and rotation period.
    Note on coördinate system: +z is north."""
    def __init__(self,mass, radius, period, atmosphere):
        self.mass=mass
        self.radius=radius
        self.period=period
        self.atmosphere=atmosphere

    def altitude(self,coordinates):
        return np.linalg.norm(coordinates)-self.radius
        
    def get_atm_rho(self,coordinates):
        return self.atmosphere.get_rho(0,self.altitude(coordinates))
        
    def gravity(self,coordinates):
        return -(coordinates/np.linalg.norm(coordinates))*(G*self.mass)/(coordinates.dot(coordinates))
        
    def absoluteWind(self,coordinates):
        """Returns velocity vector of air as it rotates. Assume that the planet
        is rotating about the z-axis. Atmosphere is assumed to be rotating as a
        solid body at the same rate as the rest of the planet up to 1000km, and 
        is then assumed to be stationary. """
        if self.altitude(coordinates)>1e6:
            return np.array([0.,0.,0.])
        else:
            r=sqrt(coordinates[0]**2+coordinates[1]**2)
            velocity = np.array([-coordinates[1],coordinates[0],0.])*((2*pi)/self.period)
            return velocity
            
    def fiveKeplerianElements(self,position,velocity):
        """Returns an array: eccentricity, semimajor axis, inclination, 
        longitude of ascending node, argument of periapsis."""
        
        r=np.linalg.norm(position)
        v=np.linalg.norm(velocity)
        
        I=np.array([1.,0.,0.])
        K=np.array([0.,0.,1.])
        
        h_vec=np.cross(position,velocity)
        h_sca=np.linalg.norm(h_vec)
        mu=G*self.mass
        
        n=np.cross(K,h_vec)
        
        a=(-mu/2)*1/(v**2/2-mu/r)
        ecc_vec=(1/mu)*(v*v*position - (position.dot(velocity))*velocity)-position/r
        ecc=np.linalg.norm(ecc_vec)
        
        inclination=angle_between(h_vec,K)
        longAscendingNode=angle_between(n,I)
        argPeriapsis=np.arccos((n.dot(ecc_vec))/(ecc*np.linalg.norm(n)))
        
        return [ecc,a,inclination,longAscendingNode,argPeriapsis]
        
    def isBelowSurface(self,position):
        r=np.linalg.norm(position)
        if r-self.radius<crash_alt:
            return True
        else:
            return False
        

class Stage(object):
    """A rocket stage.
    TODO (HIGH): more realistic I_sp vs altitude
    TODO (MEDIUM): support for C_D vs AoA and Mach rather than constant C_D
    TODO (LOW): support for nonzero L/D ratios
    """
    def __init__(self,M,M_R,vac_Isp,SL_Isp,additional_C_D,additional_frontal_area,max_throttle,min_throttle,throttle_schedule,prop_margin,tvc_schedule):
        self.drymass=M/M_R;
        self.propmass=M-self.drymass
        self.residualProp=self.propmass*prop_margin
        self.vac_Isp=vac_Isp
        self.SL_Isp=SL_Isp
        self.C_D_additional=additional_C_D
        self.additionalFrontalArea=additional_frontal_area
        self.maxThrottle=max_throttle
        self.minThrottle=min_throttle
        
        """A weird halfway thing: integrate actual propellant flow values (as 
        opposed to relative values in throttle_schedule) over unit time. This 
        essentially gives us the mean throttle value."""
        meanThrottle=min_throttle
        
        for j in range(len(throttle_schedule[0])-1):
            meanThrottle+=(max_throttle-min_throttle)*0.5*(throttle_schedule[1][j+1]+throttle_schedule[1][j])*(throttle_schedule[0][j+1]-throttle_schedule[0][j])
        
        if meanThrottle>0:
            self.burnoutTime=(self.propmass*(1-prop_margin))/meanThrottle
        else:
            self.burnoutTime=3e7
            
        self.ignitionTime=0.0
        self.sepTime=self.burnoutTime
        
        self.throttleSchedule=np.copy(throttle_schedule)
        self.tvcSchedule=np.copy(tvc_schedule)
        
        self.throttleSchedule[0][0:len(self.throttleSchedule[0])]*=self.burnoutTime
        self.throttleSchedule[1][0:len(self.throttleSchedule[1])]*=(max_throttle-min_throttle)
        self.throttleSchedule[1][0:len(self.throttleSchedule[1])]+=min_throttle
        
        self.tvcSchedule[0][0:len(self.tvcSchedule[1])]*=self.burnoutTime

        print(self.throttleSchedule)
            
        self.mdotInterpolated=spip.interp1d(self.throttleSchedule[0],self.throttleSchedule[1],kind='quadratic')
        self.TVCsInterpolated=spip.interp1d(self.tvcSchedule[0],self.tvcSchedule[1:3,:],kind='quadratic')
        
    def burnFuel(self,mass):
        self.propmass-=mass
        print("Stage.burnFuel is deprecated!")
        
    def mass(self,time):
        if time>=self.burnoutTime:
            if time>=self.sepTime:
                return 0.
            else:
                return self.drymass+self.residualProp
        elif time<self.ignitionTime:
            return self.drymass+self.propmass
        else:
            lastControlPointIndex=0
            
            for j in range(len(self.throttleSchedule[0])-1):
                if (time>=self.throttleSchedule[0][j]) and (time<self.throttleSchedule[0][j+1]):
                    lastControlPointIndex=j
                    
            propburned=0.0
            
            for j in range(lastControlPointIndex):
                propburned+=0.5*(self.throttleSchedule[1][j+1]+self.throttleSchedule[1][j])*(self.throttleSchedule[0][j+1]-self.throttleSchedule[0][j])
            
            propburned+=0.5*(self.mdotInterpolated(time)+self.throttleSchedule[1][lastControlPointIndex])*(time-self.throttleSchedule[0][lastControlPointIndex])
            
            return self.drymass+self.propmass-propburned

    def thrust(self,t,x,planet):
        """I_sp variation with altitude needs fixing."""
        if t>=self.burnoutTime or t<self.ignitionTime:
            return np.zeros(3)
        else:
            rho_SL=planet.atmosphere.rho_SL
            sl_v_e=self.SL_Isp*9.81
            vc_v_e=self.vac_Isp*9.81
            thrustMag=self.mdotInterpolated(t)*(sl_v_e+(vc_v_e-sl_v_e)*((rho_SL-planet.get_atm_rho(x))/rho_SL))
        
            thphi=self.TVCsInterpolated(t);
            unitThrustVector=np.array([sin(thphi[0])*cos(thphi[1]), sin(thphi[0])*sin(thphi[1]),cos(thphi[0])])
        
            return thrustMag*unitThrustVector
    
    def dragCoefficientContribution(self,t):
        if t<self.sepTime:
            return self.C_D_additional
        else:
            return 0
            
    def frontalAreaContribution(self,t):
        if t<self.sepTime:
            return self.additionalFrontalArea
        else:
            return 0
        
    def shiftTimesBy(self,preIgnitionCoastTime,postBurnoutCoastTime):
        self.ignitionTime+=preIgnitionCoastTime
        self.burnoutTime+=preIgnitionCoastTime
        self.sepTime+=preIgnitionCoastTime+postBurnoutCoastTime
        
        self.tvcSchedule[0][0:len(self.tvcSchedule[0])]+=preIgnitionCoastTime
        self.throttleSchedule[0][0:len(self.throttleSchedule[0])]+=preIgnitionCoastTime
        
        self.mdotInterpolated=spip.interp1d(self.throttleSchedule[0],self.throttleSchedule[1],kind='quadratic')
        self.TVCsInterpolated=spip.interp1d(self.tvcSchedule[0],self.tvcSchedule[1:3,:],kind='quadratic')        
            
class Rocket(object):
    """A rocket with an arbitrary number of stages, each of which has its own 
    mass ratio, total mass, vacuum Isp, sea-level Isp, independent coefficient 
    of drag, coefficient of drag when attached to upper stages,...
    Also, a propellant usage schedule, staging times, and a thrust vectoring 
    schedule.
        """

    def __init__(self,payload):
    #TODO (HIGH): implement fairings.
    #TODO (MEDIUM): implement parallel staging & propellant cross-feed
    #TODO (LOW): make the rocket a rocket rather than a point mass, and include rotational inertia (so tvc schedule is actual tvc schedule rather than magically pointing thrust vector elsewhere).        
        self.stages=[]
        self.flightlog=[]
        self.stages.append(copy.copy(payload))
#        self.ignitionTimes=[]
#        self.separationPropLevels=[]
#        self.propellantUsageSchedule=propellantSchedule
#        self.tvcSchedule=tvcSchedule
        
#        self.position=np.zeros(3)
#        self.velocity=np.zeros(3)
        
    def addSerialStage(self,stageToAdd,preignitionCoastTime,postBurnoutCoastTime):
        stageCopy=copy.copy(stageToAdd)
        
        stageCopy.shiftTimesBy(preignitionCoastTime,postBurnoutCoastTime)
        
        for stage in self.stages:
            stage.shiftTimesBy(stageCopy.sepTime,0.0)

        self.stages.append(stageCopy)
        
    def addStage(self,stageToAdd,separationPropLevel,ignitionTime):
        """Stage is added to "bottom" of rocket - build the rocket 
        from top down"""
        self.stages.append(copy.copy(stageToAdd))
        self.separationPropLevels.append(separationPropLevel)
        self.ignitionTimes.append(ignitionTime)
        print("deprecated!")
        
    def separateStage(self):
        self.ignitionTimes.remove(-1)
        self.separationPropLevels.remove(-1)
        print("deprecated!")
        return self.stages.pop()
        
    def totalMass(self):
        m=self.totalPropMass()
        for S in self.stages:
            m+=S.drymass
        return m
        print("deprecated!")
        
    def totalPropMass(self):
        m=0
        for S in self.stages:
            m+=S.fuelmass
        return m
        print("deprecated!")
        
    def dmdt(self,mdot,fuelmass,engineRunning):
        print("deprecated!")
        if fuelmass>0 and engineRunning:
            return min(self.stages[-1].maxThrottle,max(mdot,self.stages[-1].minThrottle))
        else:
            return 0
    
    def dvdt(self,x,v,m,thePlanet,thrust):
        """note that x,v,m may not actually be the same as those stored in self,
        because of RK4!
        returns an array of 3-arrays: a from thrust, a from gravity, a from drag
        """
        a_g=thePlanet.gravity(x)
        
        v_air=v-thePlanet.absoluteWind(x)

        airspeed=np.linalg.norm(v_air)
        a_d=-v_air*(1.0/m)*(0.5*thePlanet.get_atm_rho(x)*airspeed*self.stages[-1].C_D_below*self.stages[-1].frontalArea)

        a_t=thrust/m
        
        return [a_t, a_d, a_g]

    def v_e_atAltitude(self,x,planet):
        """This needs to be fixed later. Currently it's just about as wrong as 
        it's possible for it to be."""
        rho_SL=planet.atmosphere.rho_SL
        sl_v_e=self.stages[-1].SL_Isp*9.81
        vc_v_e=self.stages[-1].vac_Isp*9.81
        return sl_v_e+(vc_v_e-sl_v_e)*((rho_SL-planet.get_atm_rho(x))/rho_SL)
    
    def sDot(self,t,s,planet,targetSMA):
        """s is the rocket's state vector, consisting of: x_{1,2,3},v_{1,2,3},m.
        Note that a and m are purely a function of t."""

        x=s[0:3]
        v=s[3:6]
        
        kepEls=planet.fiveKeplerianElements(x,v)
        
        sdot=np.zeros(7)
        
        if kepEls[1]>=targetSMA:
            sdot[6]=0
            #print(s)
            print(str(kepEls[1])+">="+str(targetSMA))
            return sdot            
        else:
            sdot[0:3]=v
            
            F_t=self.thrust(t,x,planet)
            a_t=F_t/self.mass(t)
            
            v_air=v-planet.absoluteWind(x)
            airspeed=np.linalg.norm(v_air)
            a_d=-v_air*(0.5*planet.get_atm_rho(x)*airspeed*self.dragCoefficient(t)*self.frontalArea(t))/self.mass(t)
            
            a_g=planet.gravity(x)
            
            sdot[3:6]=a_t+a_d+a_g
            
            sdot[6]=1
            
            flightstatus=np.zeros(5)
            flightstatus[0]=t
            flightstatus[1]=self.mass(t)
            flightstatus[2]=np.linalg.norm(a_t)
            flightstatus[3]=np.linalg.norm(a_d)
            flightstatus[4]=np.linalg.norm(a_g)
                                    
            self.flightlog.append(flightstatus)
                        
            return sdot
            

    def thrust(self,t,x,planet):
        """Sum thrusts of all stages at t. This allows for future features like parallel staging."""
        thrust=np.zeros(3)
        for stage in self.stages:
            thrust+=stage.thrust(t,x,planet)
        return thrust
        
    def mass(self,t):
        mass=0.
        for stage in self.stages:
            mass+=stage.mass(t)
        return mass
            
    def dragCoefficient(self,t):
        drag_coefficient=0.
        for stage in self.stages:
            drag_coefficient+=stage.dragCoefficientContribution(t)
        return drag_coefficient
        
    def frontalArea(self,t):
        frontal_area=0.
        for stage in self.stages:
            frontal_area+=stage.frontalAreaContribution(t)
        return frontal_area
    
    def integrateFlight(self,initialCoords,initialVel,planet,t0,maxSimTime, stopOnSuccess, desiredSemiMajor):
        #TODO: add support for stopOnSuccess
        backend='dopri5'
        
        s0=np.zeros(7)
        s0[0:3]=initialCoords
        s0[3:6]=initialVel
        s0[6]=1
        
        #solver = odeint(self.sDot).set_integrator(backend)
        #solver.set_initial_value(s0, t0)
        #solver.set_f_params(planet,desiredSemiMajor)
        
        return odeint(self.sDot,s0,t0,args=(planet,desiredSemiMajor),full_output=True)
        
    def rk4Flight(self,initialCoords,initialVel,planet,t0,simTime,stopOnSuccess,desiredSemiMajor,h_base):

        singularitytimes=[]
        for s in self.stages:
            singularitytimes.append(s.ignitionTime)
            singularitytimes.append(s.burnoutTime)
            singularitytimes.append(s.sepTime)
        
        h_min=1e-12
        
        terminalGuidanceThreshold=200000*h_base
        singularityThreshold=3*h_base
        terminalGuidanceActive=False
        
        self.flightlog=[]

        s=np.zeros(7)
        s[0:3]=initialCoords
        s[3:6]=initialVel
        s[6]=1
        
        t=t0
        
        traj=[]
        
        keepSimulating=1
        
        while t<simTime and keepSimulating>0 and t<self.stages[1].burnoutTime:
            distanceToNextSingularity=3e7
            
            for time in singularitytimes:
                if t<time and time-t<distanceToNextSingularity:
                    distanceToNextSingularity=time-t
            
            traj.append(np.copy(s))
            k1=self.sDot(t,s,planet,desiredSemiMajor)

            kepEls=planet.fiveKeplerianElements(s[0:3],s[3:6])
            
            if desiredSemiMajor-kepEls[1]<terminalGuidanceThreshold:
                h=max(h_min, h_base*(desiredSemiMajor-kepEls[1])/terminalGuidanceThreshold)
                if not terminalGuidanceActive:
                    terminalGuidanceActive=True
                    print("Terminal guidance active!")
            elif distanceToNextSingularity<singularityThreshold:
                h=max(h_min, h_base*distanceToNextSingularity/singularityThreshold)
            else:
                h=h_base

            k2=self.sDot(t+0.5*h,s+(0.5*h)*k1,planet,desiredSemiMajor)
            k3=self.sDot(t+0.5*h,s+(0.5*h)*k2,planet,desiredSemiMajor)            
            k4=self.sDot(t+h,s+h*k3,planet,desiredSemiMajor)
            
            s+=h*(k1+2*k2+2*k3+k4)/6
            t+=h
            
            if(stopOnSuccess==True):
                keepSimulating=k4[6]
         
        return traj
         
    def simulateFlight(self,initialCoords,initialVel, planet, maxSimTime, timeStep, stopAtLastMeco, desiredSemiMajor):
        """Simulate the rocket's flight. Assumes all stages have been added and 
        that the propellant and TVC schedules have been defined. Also assumes
        that initialVel is initial velocity relative to ground/atmosphere.
            stopAtLastMeco: terminate simulation at final stage MECO.
        Returns a list of arrays containing the following rows:
            time
            position (three)
            velocity (three)
            acceleration from thrust (three)
            acceleration from drag (three)
            acceleration from gravity (three)
            Q
            remaining propellant"""
        
        x=np.copy(initialCoords)
        v=np.copy(initialVel)
        v+=planet.absoluteWind(x)
        
        interpPropSchedule=spip.interp1d(self.propellantUsageSchedule[0],self.propellantUsageSchedule[1],kind='linear')
        interpTVCSchedule=spip.interp1d(self.tvcSchedule[0],self.tvcSchedule[1:3,:],kind='linear')
        
        t=0.0
        eventTimer=0.0
        nextEventType=0 #0 is engine ignition, 1 is ECO and stage sep
        
        enginesRunning=True
        stopSimulation=False
        
        allStates=[]
        thphi=interpTVCSchedule(t);
        unitThrustVector=np.array([sin(thphi[0])*cos(thphi[1]), sin(thphi[0])*sin(thphi[1]),cos(thphi[0])])
        thrust=self.v_e_atAltitude(initialCoords,planet)*interpPropSchedule(t)*unitThrustVector
        startAcceleration=self.dvdt(x,v,self.totalMass(),planet,thrust)
        
        """drawn out this way for improved execution speed. I'm sure there's a 
        better way to do this..."""
        vehState=np.zeros(19)
        vehState[0]=t
        vehState[1:4]=np.copy(x)
        vehState[4:7]=np.copy(v)
        vehState[7:10]=startAcceleration[0]  #a_thrust
        vehState[10:13]=startAcceleration[1] #a_drag
        vehState[13:16]=startAcceleration[2] #a_gravity
        vehState[16]=0.0
        vehState[17]=self.totalPropMass()
        vehState[18]=self.totalMass()
        
        allStates.append(vehState)
              
        while t<maxSimTime and stopSimulation==False:
            #print("---DEBUG: X AND V---")
            #print(x)
            #print(v)
            
            #Engine and staging logic
            if nextEventType==0 and eventTimer>=self.ignitionTimes[-1]:
                self.ignitionTimes.pop()
                enginesRunning=True
                nextEventType=1
                eventTimer=0.0
                print("ignition, ", end="")
                
            if nextEventType==1 and self.stages[-1].fuelmass<=self.separationPropLevels[-1] and len(self.stages)>1:
                self.separationPropLevels.pop()
                discardedStage=self.stages.pop()
                enginesRunning=False
                nextEventType=0
                eventTimer=0.0
                print("MECO & stage sep (t="+str(t)+", alt="+str(planet.altitude(x))+") ", end="")
            
            if len(self.stages)==1 and self.stages[-1].fuelmass<=self.separationPropLevels[-1]:
                enginesRunning=False
                print("upper stage burnout! t="+str(t)+", alt="+str(planet.altitude(x)))
                if stopAtLastMeco:
                    stopSimulation=True
            
            #decide on dt
            thphi=interpTVCSchedule(t);
            unitThrustVector=np.array([sin(thphi[0])*cos(thphi[1]), sin(thphi[0])*sin(thphi[1]),cos(thphi[0])])
            thrust=self.v_e_atAltitude(initialCoords,planet)*interpPropSchedule(t)*unitThrustVector            
            a_guess=self.dvdt(x,v,self.totalMass(),planet,thrust)
            a_mag_guess=np.linalg.norm(a_guess[0])+np.linalg.norm(a_guess[1])+np.linalg.norm(a_guess[2])
            
            dt=timeStep/a_mag_guess
            
            #RK4, using for-loop this time. 
            km=np.zeros(4)
            kv=np.zeros([4,3])
            kx=np.zeros([4,3])

            #print("---START RK---")
            for j in range(4):
                """It's nice that Python lets you use -1 as an index. Makes this
                bit a lot cleaner."""
                substep=0.5*floor((j+1)/2)
                
                km[j]=self.dmdt(interpPropSchedule(t+substep*dt),self.stages[-1].fuelmass-dt*substep*km[j-1],enginesRunning)
                
                thphi=interpTVCSchedule(t+substep*dt);
                unitThrustVector=np.array([sin(thphi[0])*cos(thphi[1]), sin(thphi[0])*sin(thphi[1]),cos(thphi[0])])
                thrust=self.v_e_atAltitude(x,planet)*km[j]*unitThrustVector
                a=self.dvdt(x+substep*dt*kx[j-1,:], v+substep*dt*kv[j-1,:], self.totalMass()-substep*dt*km[j-1],planet,thrust)
                
                #print(a)
                
                kv[j,:]=a[0]+a[1]+a[2]
                
                kx[j,:]=v+substep*dt*kv[j-1]  

            #print("---RK KONSTANTS---")            
            
            #print(km)   
            #print(kv)
            #print(kx)
            
            #print(dt)
            
            x+=dt*(kx[0]+2*kx[1]+2*kx[2]+kx[3])/6
            v+=dt*(kv[0]+2*kv[1]+2*kv[2]+kv[3])/6
            
            self.stages[-1].burnFuel(dt*(km[0]+2*km[1]+2*km[2]+km[3])/6)
            
            #stop if we've crashed
            if planet.isBelowSurface(x):
                stopSimulation=True
                print('crashed!')
            
            #stop if we have enough orbital energy
            semiMajor=planet.fiveKeplerianElements(x,v)[1]
            if semiMajor>=desiredSemiMajor and enginesRunning:
                enginesRunning=False
                print("Enough energy! t="+str(t)+", alt="+str(planet.altitude(x)))
                if stopAtLastMeco:
                    stopSimulation=True

            t+=dt
            eventTimer+=dt
                        
            #Do appending and stuff
            vehState=np.zeros(19)
            vehState[0]=t
            vehState[1:4]=x
            vehState[4:7]=v
            vehState[7:10]=a[0]
            vehState[10:13]=a[1]
            vehState[13:16]=a[2]
            vehState[16]=0.5*(np.linalg.norm(v-planet.absoluteWind(x))**2)*planet.get_atm_rho(x)
            vehState[17]=self.totalPropMass()
            vehState[18]=self.totalMass()
            
            allStates.append(vehState)
            
        
        return allStates
        
class Simulator(object):
    def __init__(self,planet,rocket,propellantSchedule,tvcSchedule):
        self.tvcSchedule=tvcSchedule
        self.propellantSchedule=propellantSchedule
        self.rocket=rocket
        self.planet=planet
    
    def plfitness(self,initialCoords,initialVel,maxSimTime, timeStep, stopAtLastMeco, desiredSemiMajor,desiredEccentricity):
        flightStats=self.rocket.simulateFlight(initialCoords,initialVel, self.planet, maxSimTime, timeStep, stopAtLastMeco, desiredSemiMajor);
        arrayFlightStats=np.asarray(flightStats).T
        
        eccTol=0.002
        smaTol=15000
        
        smaWeight=0.001
        eccWeight=20
        propResidualWeight=-0.1
        
        kepEls=self.planet.fiveKeplerianElements(arrayFlightStats[1:4,-1],arrayFlightStats[4:7,-1])
        propResidual=arrayFlightStats[17, -1]-self.rocket.separationPropLevels[-1]
        
        fitness=0.;
        componentstr=""
        
        if abs(desiredSemiMajor-kepEls[1])>smaTol:
            fitness+=(abs(desiredSemiMajor-kepEls[1])-smaTol)*smaWeight
            componentstr+="a was "+str(kepEls[1])+", "
            
        if abs(kepEls[0]-desiredEccentricity)>eccTol:
            fitness+=(abs(kepEls[0]-desiredEccentricity)-eccTol)*eccWeight
            componentstr+="e was "+str(kepEls[0])
        else:
            fitness+=propResidual*propResidualWeight
            componentstr+="e good and "+str(propResidualWeight)+"kg residual."
        
        print("error="+str(fitness)+" ("+componentstr+")")
        return fitness;