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

testStartPoint=np.array([1.0,1.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.5,0,0.4,0.8,0.5,0.8,1.2,1.6,0.5,0.5,0.5,0.5,0.5,0.5,1.0])

#testStartPoint=np.asarray([ 0.97368009,  0.96763239,  0.24523624,  0.49324878,  0.72762755, 1.00416857,  0.95933795,  0.97403668,  0.95350731,  1.00531046, 0.46458082,  1.67220112,  2.17312733,  2.66346563,  0.41540276, 2.4781154 ,  3.13789563,  3.55330601])
#testStartPoint=np.array([  3.11592812,   2.31278171,   0.25      ,   0.5       , 0.75      ,   3.98778947,   4.19956555,   1.90349192, 3.27566526,   2.63714774, -14.90366124,   4.38034555, 21.22269636,  23.06061613,  -1.72499083,   7.37530809, 5.48134183,  -2.59405555])
#testStartPoint=np.array([  1,   1,   0.25      ,   0.5       , 0.75      ,   1,   1,   1, 1,   1, -14.90366124,   4.38034555, 21.22269636,  23.06061613,  -1.72499083,   7.37530809, 5.48134183,  -2.59405555])

def f(pars):
    """for the current test, we have 0+2 upperstage throttle, 3+4 lowerstage throttle,
    1+3 upperstage tvc-theta and 1+3 lowerstage tvc-theta points - 17 in all"""
    
    upperThrottleTimes=[0.,0.5,1.]
    upperThrottleSettings=[np.clip(pars[0],0.,1.),np.clip(pars[23],0.,1.),np.clip(pars[1],0.,1.)]
    
    lowerThrottleTimes=[0.,0.25+np.clip(pars[2],-0.125,0.125),0.5+np.clip(pars[3],-0.125,0.125),0.75+np.clip(pars[4],-0.125,0.125),1.]
    lowerThrottleSettings=[1.,np.clip(pars[5],0.,1.),np.clip(pars[6],0.,1.),np.clip(pars[7],0.,1.),np.clip(pars[8],0.,1.)]    
    
    upperThrottleSchedule=np.array([upperThrottleTimes,upperThrottleSettings])
    lowerThrottleSchedule=np.array([lowerThrottleTimes,lowerThrottleSettings])
    
    lowerTVCPhi=[np.clip(pars[10],0.,pi)+pi/2,np.clip(pars[11],0.,pi)+pi/2,np.clip(pars[11],0.,pi)+pi/2,np.clip(pars[11],0.,pi)+pi/2,np.clip(pars[12],0.,pi)+pi/2]
    upperTVCPhi=[np.clip(pars[14],0.,pi)+pi/2,np.clip(pars[20],0.,pi)+pi/2,np.clip(pars[21],0.,pi)+pi/2,np.clip(pars[15],0.,pi)+pi/2,np.clip(pars[16],0.,pi)+pi/2]
    

    lowerTVCSchedule=np.array([[0.,np.clip(pars[ 9],0.,1.),np.clip(pars[16],0.,1.),np.clip(pars[17],0.,1.),1.0],[pi/2,pi/2,pi/2,pi/2,pi/2],lowerTVCPhi])
    upperTVCSchedule=np.array([[0.,np.clip(pars[13],0.,0.33),np.clip(pars[18],0.33,0.66),np.clip(pars[19],0.66,1.),1.],[pi/2,pi/2,pi/2,pi/2,pi/2],upperTVCPhi])

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

test=spop.fmin(f,testStartPoint,maxiter=10)
#last time: test=np.asarray([  9.18614021e-01,   8.91426376e-01,   1.70268538e-04, -1.06837974e-03,  -2.04178832e-04,   9.59172540e-01, 8.67570250e-01,   1.02605773e+00,   1.02985045e+00, 3.55660691e-01,  -4.23176627e-04,   4.16105945e-01, 1.00786586e+00,   4.98887445e-01,   8.54004016e-01, 1.67364199e+00,   2.06206756e+00])

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