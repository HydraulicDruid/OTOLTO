from otolto import * 
import pylab

"""Information that could be useful to plot:
    Velocity
    Total acceleration
    Orbital diagram
    Q
    Altitude
    AoA
Information that could be useful to calculate:
    Delta-V loss to gravity drag
    Delta-V loss to airdrag
Hillclimbing information to display:
    Total iterations, missteps made in this cycle
    % of missteps
    Fitness vs last best fitness
    Excess propellant, orbital elements
"""

def rnd():
    return 2*(np.random.rand()-0.5)

eccentricityWeightCap=10000.0
inclinationWeightCap=0.1
excessFuelWeightingFactor=0.0 #per kilogram - 200 is unusually high, of course

desiredEccentricity=0.0
desiredInclination=0.0
desiredSemiMajor=7e6

maximumAttempts=20

testAtmosphere=Atmosphere(0,1.2,29.26,260)
testEarth=Planet(5.97219e24, 6371e3, 86164.1, testAtmosphere)


testUpperStage=Stage(120.0,7.0,350.0,0.0,0.7,0.7,0.28,0.3,0.3)
testLowerStage=Stage(1400.0,7.2,290.0,265.0,0.7,0.8,0.28,7.0,4.0)

propellantSchedule=np.array([[0.0,175.0,360.0,3e7],[7.0,5.0,0.3,0.3]])

#rows are time, theta, phi (theta is azimuth and phi is polar)
tvcSchedule=np.array([[0.,15.,50.,200.,600.,3e7], [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2],                                     [pi/2,pi/2,pi/2,pi,pi,pi]])

misstepsBeforeReduction=10
reductionStep=2.0
reductionFactor=1.0
misstepsMade=0
oldFitness=0.0

esT=np.linspace(-pi,pi,2000)
esX,esY=6371e3*np.sin(esT),6371e3*np.cos(esT)
pylab.plot(esX,esY)

minX=1e9
maxX=-1e9
minY=1e9
maxY=-1e9

for j in range(maximumAttempts):
    #A horrible kludge, of course. Needs fixing.
    propellantScheduleRand=np.array([[0.0,120*rnd(),120*rnd(),0.0],[3*rnd(),3*rnd(),0.1*rnd(),0.1*rnd()]])
    tvcScheduleRand=np.array([[0.,5*rnd(),25*rnd(),100*rnd(),300*rnd(),0.],[0.,0.,0.,0.,0.,0.],[0.0,0.8*rnd(),0.8*rnd(),0.8*rnd(),0.8*rnd(),0.0]])

    oldPropellantSchedule=np.copy(propellantSchedule)
    oldTVCSchedule=np.copy(tvcSchedule)
    
    propellantSchedule+=propellantScheduleRand*reductionFactor
    tvcSchedule+=tvcScheduleRand*reductionFactor
    
    testRocket=Rocket(propellantSchedule,tvcSchedule)
    testRocket.addStage(testUpperStage,1,5)
    testRocket.addStage(testLowerStage,20,0)
 
    coords=np.array([0.0, 6371010.0, 0.0])
    velocity=np.array([0.0,0.0,0.0])
    #print(str(coords)+str(velocity))

    #yes, I am extremely professional and name my variables "blah". (TODO: change this)
    blah=testRocket.simulateFlight(coords,velocity,testEarth,600.0,8.0,True,desiredSemiMajor)
    arblah=np.asarray(blah).T

    minX=min(min(arblah[1]),minX)
    maxX=max(max(arblah[1]),maxX)    

    minY=min(min(arblah[2]),minY)
    maxY=max(max(arblah[2]),maxY)

    xwidth=maxX-minX
    yheight=maxY-minY

    pylab.xlim( minX-0.1*xwidth, maxX+0.1*xwidth )
    pylab.ylim( minY-0.1*yheight, maxY+0.1*yheight )
    pylab.plot(arblah[1],arblah[2])

    kepEls=testEarth.fiveKeplerianElements(arblah[1:4,-1],arblah[4:7,-1])

    #print(arblah[17,-1])

    fitness=arblah[17, -1]*excessFuelWeightingFactor + min(eccentricityWeightCap, 1/abs(kepEls[0]-desiredEccentricity)) +   min(inclinationWeightCap, 1/abs(kepEls[2]-desiredInclination))

    print("iteration! fitness="+str(fitness)+" from weighting factors "+str(kepEls[0])+","+str(kepEls[2])+" degrees, "+str(arblah[17, -1])+"kg.")
    if fitness>oldFitness:
        misstepsMade=0
        oldFitness=fitness
        print("better trajectory!")
    else:
        propellantSchedule=np.copy(oldPropellantSchedule)
        tvcSchedule=np.copy(oldTVCSchedule)
        misstepsMade+=1
        
    if misstepsMade>=misstepsBeforeReduction:
        reductionFactor/=reductionStep
        print('reduced step size!')
    print("---------------------------------")    
    
    del(testRocket)
    del(propellantScheduleRand)
    del(tvcScheduleRand)
    del(coords)
    del(velocity)