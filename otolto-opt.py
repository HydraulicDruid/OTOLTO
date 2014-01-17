from otolto import * 
import scipy.optimize as spop

eccentricityWeightCap=10000.0
inclinationWeightCap=1.0
excessFuelWeightingFactor=0.0 #per kilogram - 200 is unusually high, of course

desiredEccentricity=0.0
desiredInclination=0.0
desiredSemiMajor=7e6

testAtmosphere=Atmosphere(0,1.2,29.26,260)
testEarth=Planet(5.97219e24, 6371e3, 86164.1, testAtmosphere)

testUpperStage=Stage(120.0,7.0,350.0,0.0,0.7,0.7,0.28,0.4,0.2)
testLowerStage=Stage(1400.0,7.2,290.0,265.0,0.7,0.8,0.28,7.0,4.0)

coords=np.array([0.0, 6371010.0, 0.0])
velocity=np.array([0.0,0.0,0.0])

steps=10

arr=np.zeros((steps,steps))

def f(blah):
    propellantSchedule=np.array([[0.0,3e7],[7.667,0.405]])

    #ay=2*(j/(steps-1))[
    #be=2*(k/(steps-1))
    #print(str(ay)+","+str(be))
    #rows are time, theta, phi (theta is azimuth and phi is polar)
    tvcSchedule=np.array([[0.,blah[3]*100.,blah[4]*100.,blah[5]*100.,3e7], [pi/2,pi/2,pi/2,pi/2,pi/2],[pi/2,pi/2+blah[0],pi/2+blah[1],pi/2+blah[2],pi/2+blah[2]]])

    testRocket=Rocket(propellantSchedule,tvcSchedule)
    testRocket.addStage(testUpperStage,1,5)
    testRocket.addStage(testLowerStage,20,0)
        
    testSim=Simulator(testEarth,testRocket,propellantSchedule,tvcSchedule)
        
    #arr[j][k]=testSim.plfitness(coords,velocity,700,10., True, desiredSemiMajor)

    return testSim.plfitness(coords,velocity,700,10., True, desiredSemiMajor, desiredEccentricity)
        
    #del(testRocket)
    #del(testSim)
        
    #del(propellantSchedule)
    #del(tvcSchedule)

initialThing=np.zeros((2,3))
for k in range(3):
    initialThing[0][k]=np.random.rand()*0.8
    initialThing[1][k]=np.random.rand()*2.0
    
initialThing=np.cumsum(initialThing, axis=1)
initialThing=np.concatenate((initialThing[0],initialThing[1]), axis=1)

print(initialThing)

#spop.fmin_bfgs(f, [-0.1076601, 2.10792968, 0.37004736, 4.03], maxiter=2)
test=spop.fmin(f, [  2.44515213e-04,   9.67057808e-01,   2.10612178e+00, 3.31306112e-01,   2.24395094e+00,   3.98517198e+00], maxiter=10)
#print("closest found: "+str(test))
    
#for j in range(steps):
#    for k in range(steps):
#        ay=-0.15+0.3*(j/(steps-1))
#        be=1.9+0.3*(k/(steps-1))
#        print(str(ay)+","+str(be))
#        arr[j,k]=f([ay,be,0.37,4.])