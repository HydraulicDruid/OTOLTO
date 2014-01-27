from otolto import * 
import scipy.optimize as spop
import cma
import multiprocessing

if __name__ == '__main__':
    print('My name is __main__!')
    multiprocessing.freeze_support()

multiprocessing.freeze_support()

desiredEccentricity=0.0
desiredInclination=0.0
desiredSemiMajor=7e6

eccTol=0.002
smaTol=1.0
        
smaWeight=1e-4
eccWeight=20
propResidualWeight=-10.

testAtmosphere=Atmosphere(0,1.2,29.26,260)
testEarth=Planet(5.97219e24, 6371e3, 86164.1, testAtmosphere)

coords=np.array([0.0, 6371010.0, 0.0])
velocity=np.array([0.0,0.0,0.0])

steps=10

arr=np.zeros((steps,steps))

testSPoint=np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.0,0.5,0.1,0.2,0.7,0.4,0.1,0.2,.35,7.0,0.5])
testLBound=np.array([.01,.01,.01,.01,.01,.01,.01,.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.5,0.0])
testUBound=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, pi, pi, pi, pi, pi, pi, pi, pi,.50,11.,99.])

#testStartPoint=np.asarray([ 0.97368009,  0.96763239,  0.24523624,  0.49324878,  0.72762755, 1.00416857,  0.95933795,  0.97403668,  0.95350731,  1.00531046, 0.46458082,  1.67220112,  2.17312733,  2.66346563,  0.41540276, 2.4781154 ,  3.13789563,  3.55330601])
#testStartPoint=np.array([  3.11592812,   2.31278171,   0.25      ,   0.5       , 0.75      ,   3.98778947,   4.19956555,   1.90349192, 3.27566526,   2.63714774, -14.90366124,   4.38034555, 21.22269636,  23.06061613,  -1.72499083,   7.37530809, 5.48134183,  -2.59405555])
#testStartPoint=np.array([  1,   1,   0.25      ,   0.5       , 0.75      ,   1,   1,   1, 1,   1, -14.90366124,   4.38034555, 21.22269636,  23.06061613,  -1.72499083,   7.37530809, 5.48134183,  -2.59405555])

def f_TVCOnly(pars):
    """pars has 4n+3 datapoints (where n is the number of data points, and the last two points are the throttles of the stages, and the actual last one is the upper stage coast time)"""
    
    #upperThrottleTimes=[0.,1.]
    #upperThrottleSettings=[np.clip(pars[0],0.,1.),np.clip(pars[1],0.,1.)]
    nonDataPointCount=3
    
    throttleTimes=[0.,1.]
    throttleSettings=[1.,1.]
    
    upperThrottleSchedule=np.array([throttleTimes,throttleSettings])
    lowerThrottleSchedule=np.array([throttleTimes,throttleSettings])
    
    n=floor((len(pars)-nonDataPointCount)/4)
    #print(n)
    
    lowerTVCTimes=pars[0:n]
    upperTVCTimes=pars[n:2*n]
    
    lowerTVCPhi=pars[2*n:3*n]
    upperTVCPhi=pars[3*n:4*n]
    
    lowerTVCTheta=np.ones(n+1)*(pi/2)
    upperTVCTheta=np.ones(n+1)*(pi/2)
    
    """print(lowerTVCTimes)
    print(upperTVCTimes)
    print(lowerTVCPhi)
    print(upperTVCPhi)"""
    
    #now process them
    lowerTVCTimes=np.insert(np.cumsum(np.clip(lowerTVCTimes,0.01,1)),0,0)
    lowerTVCTimes/=lowerTVCTimes[-1]
    
    upperTVCTimes=np.insert(np.cumsum(np.clip(upperTVCTimes,0.01,1)),0,0)
    upperTVCTimes/=upperTVCTimes[-1]
    
    lowerTVCPhi=np.clip(np.cumsum(np.clip(np.insert(lowerTVCPhi,0,0),0.,np.Inf)),0,pi)+pi/2
    upperTVCPhi=np.clip(np.cumsum(np.clip(np.insert(upperTVCPhi,0,lowerTVCPhi[-1]-pi/2),0.,np.Inf)),0,pi)+pi/2

    #and assemble them
    lowerTVCSchedule=np.array([lowerTVCTimes,lowerTVCTheta,lowerTVCPhi])
    upperTVCSchedule=np.array([upperTVCTimes,upperTVCTheta,upperTVCPhi])

    #print(lowerTVCSchedule)
    #print(upperTVCSchedule)

    testPayload=InertStage(5.,0.,0.,30.0)
    testUpperStage=Stage( 90.0, 8.0,350.0,  0.0,0.7,0.3,pars[-3],pars[-3],upperThrottleSchedule,0.04,upperTVCSchedule)
    testLowerStage=Stage(800.0,15.0,315.0,275.0,0.1,0.0,pars[-2],pars[-2],lowerThrottleSchedule,0.04,lowerTVCSchedule)

    testRocket=Rocket(testPayload)
    testRocket.addSerialStage(testUpperStage,np.clip(pars[-1],0.1,np.Inf),2.0)
    testRocket.addSerialStage(testLowerStage,0.0,2.0)
    
    #print("Simulating flight...")
    blah=testRocket.rk4Flight(coords,velocity,testEarth,0.0,800,True,desiredSemiMajor,3)
    arblah=np.asarray(blah)
    arlog=np.asarray(testRocket.flightlog)
    #print("done!")

    kepEls=testEarth.fiveKeplerianElements(arblah[-1,0:3],arblah[-1,3:6])
    propResidual=arlog[-1][1]-testRocket.stages[0].drymass-testRocket.stages[1].drymass-testRocket.stages[1].residualProp
    
    #print(kepEls)
    #print(propResidual)

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
        componentstr+="e good and "+str(propResidual)+"kg residual."    
    
    componentstr+=" Fitness="+str(fitness)
        
    print(componentstr)
    return fitness

"""
es=cma.CMAEvolutionStrategy(testStartPoint,1, {'maxiter':400, 'bounds':[testLBound,testUBound]})
pool=multiprocessing.Pool(es.popsize)
while not es.stop():
    X=es.ask()
    MapRes=pool.map_async(f_TVCOnly, X)
    es.tell(X, MapRes.get())
    es.disp()
"""

#test=cma.fmin(f_TVCOnly, testStartPoint, 1, maxiter=20, bounds=[testLBound,testUBound])

"""test=spop.fmin(f,testStartPoint,maxiter=100)
#last time: test=np.array([  1.62429868e-01,   4.05661539e-01,   1.70682155e-01,
         4.26791700e-01,   3.97759810e-01,   4.12189600e-01,
         6.41961575e-01,   4.37625495e-01,   2.06544490e-03,
         6.55420902e-01,   1.02260684e-01,   2.66786277e-01,
         2.26000646e-04,   8.08809548e-01,   1.22724676e-01,
         2.90182950e-01,   2.26687414e-01,   8.77939372e+00,
         7.76663207e-01])

tempPoint=np.copy(test)

nearbyFitnesses=[]

for k in range(20):
    for n in testStartPoint:
        pass
    for n in range(len(testStartPoint)):
        tempPoint[n]=test[n]+np.random.normal(scale=0.05)
    nearbyFitnesses.append(f(tempPoint))"""

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