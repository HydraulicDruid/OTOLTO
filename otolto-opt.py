from otolto import * 
import scipy.optimize as spop

desiredEccentricity=0.0
desiredInclination=0.0
desiredSemiMajor=7e6

eccTol=0.002
smaTol=1.0
        
smaWeight=0.001
eccWeight=20
propResidualWeight=-0.1

testAtmosphere=Atmosphere(0,1.2,29.26,260)
testEarth=Planet(5.97219e24, 6371e3, 86164.1, testAtmosphere)

coords=np.array([0.0, 6371010.0, 0.0])
velocity=np.array([0.0,0.0,0.0])

steps=10

arr=np.zeros((steps,steps))

testStartPoint=np.array([1.0,1.0,0.25,0.5,0.75,1.0,1.0,1.0,1.0,1.0,0.5,pi/2,pi/2+0.4,pi/2+0.8,0.5,pi/2+0.8,pi/2+1.2,pi])

testStartPoint=np.asarray([ 0.97368009,  0.96763239,  0.24523624,  0.49324878,  0.72762755, 1.00416857,  0.95933795,  0.97403668,  0.95350731,  1.00531046, 0.46458082,  1.67220112,  2.17312733,  2.66346563,  0.41540276, 2.4781154 ,  3.13789563,  3.55330601])
#testStartPoint=np.array([  3.11592812,   2.31278171,   0.25      ,   0.5       , 0.75      ,   3.98778947,   4.19956555,   1.90349192, 3.27566526,   2.63714774, -14.90366124,   4.38034555, 21.22269636,  23.06061613,  -1.72499083,   7.37530809, 5.48134183,  -2.59405555])
#testStartPoint=np.array([  1,   1,   0.25      ,   0.5       , 0.75      ,   1,   1,   1, 1,   1, -14.90366124,   4.38034555, 21.22269636,  23.06061613,  -1.72499083,   7.37530809, 5.48134183,  -2.59405555])

def f(pars):
    """for the current test, we have 0+2 upperstage throttle, 3+5 lowerstage throttle,
    1+3 upperstage tvc-theta and 1+3 lowerstage tvc-theta points - 18 in all"""
    upperThrottleSchedule=np.array([[0.,1.],[pars[0],pars[1]]])
    lowerThrottleSchedule=np.array([[0.,pars[2],pars[3],pars[4],1.],[pars[5],pars[6],pars[7],pars[8],pars[9]]])

    lowerTVCSchedule=np.array([[0.,pars[10],1.],[pi/2,pi/2,pi/2],[pars[11],pars[12],pars[13]]])
    upperTVCSchedule=np.array([[0.,pars[14],1.],[pi/2,pi/2,pi/2],[pars[15],pars[16],pars[17]]])

    testPayload=Stage(5.,1.,0.,0.,0.,0.,0.,0.,upperThrottleSchedule,0.,upperTVCSchedule)
    testUpperStage=Stage(120.0,7.0,350.0,0.0,0.7,0.28,0.4,0.2,upperThrottleSchedule,0.04,upperTVCSchedule)
    testLowerStage=Stage(1400.0,10.0,315.0,275.0,0.1,0.0,7.0,4.0,lowerThrottleSchedule,0.04,lowerTVCSchedule)

    testRocket=Rocket(testPayload)
    testRocket.addSerialStage(testUpperStage,0.5,2.0)
    testRocket.addSerialStage(testLowerStage,0.0,2.0)
        
    print("Simulating flight...")
    blah=testRocket.rk4Flight(coords,velocity,testEarth,0.0,800,True,desiredSemiMajor,0.8)
    arblah=np.asarray(blah)
    arlog=np.asarray(testRocket.flightlog)
    print("done!")

    kepEls=testEarth.fiveKeplerianElements(arblah[-1,0:3],arblah[-1,3:6])
    propResidual=arlog[-1]-testRocket.stages[0].drymass-testRocket.stages[1].drymass
    
    print(kepEls)
    print(propResidual)

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
    
    componentstr+=" Fitness="+str(fitness)
        
    print(componentstr)
    return fitness

test=spop.fmin(f,testStartPoint,maxiter=100)
#last time: test=np.asarray([ 0.97368009,  0.96763239,  0.24523624,  0.49324878,  0.72762755, 1.00416857,  0.95933795,  0.97403668,  0.95350731,  1.00531046, 0.46458082,  1.67220112,  2.17312733,  2.66346563,  0.41540276, 2.4781154 ,  3.13789563,  3.55330601])

tempPoint=np.copy(test)

nearbyFitnesses=[]

for k in range(20):
    for n in testStartPoint:
        pass
    for n in range(len(testStartPoint)):
        tempPoint[n]=test[n]+np.random.normal(scale=0.05)
    nearbyFitnesses.append(f(tempPoint))

#initialThing=np.zeros((2,3))
#for k in range(3):
#    initialThing[0][k]=np.random.rand()*0.8
#    initialThing[1][k]=np.random.rand()*2.0
    
#initialThing=np.cumsum(initialThing, axis=1)
#initialThing=np.concatenate((initialThing[0],initialThing[1]), axis=1)

#print(initialThing)

#spop.fmin_bfgs(f, [-0.1076601, 2.10792968, 0.37004736, 4.03], maxiter=2)
#test=spop.fmin(f, [  2.44515213e-04,   9.67057808e-01,   2.10612178e+00, 3.31306112e-01,   2.24395094e+00,   3.98517198e+00], maxiter=10)
#print("closest found: "+str(test))
    
#for j in range(steps):
#    for k in range(steps):
#        ay=-0.15+0.3*(j/(steps-1))
#        be=1.9+0.3*(k/(steps-1))
#        print(str(ay)+","+str(be))
#        arr[j,k]=f([ay,be,0.37,4.])