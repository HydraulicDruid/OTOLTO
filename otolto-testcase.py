from otolto import * 
import matplotlib.pyplot as plt

testAtmosphere=Atmosphere(0,1.2,29.26,260)
testEarth=Planet(5.97219e24, 6371e3, 86164.1, testAtmosphere)
testTVC=np.array([[0.,1.],[pi/2,pi/2],[pi/2+0.8,pi]])
testThrottle=np.array([[0.,1.],[1.0,1.0]])
testTVCLower=np.array([[0.,1.],[pi/2,pi/2],[pi/2,pi/2+0.8]])

testPayload=Stage(5.,1.,0.,0.,0.,0.,0.,0.,testThrottle,0.,testTVC)
testUpperStage=Stage(120.0,7.0,350.0,0.0,0.7,0.28,0.4,0.2,testThrottle,0.04,testTVC)
testLowerStage=Stage(1400.0,10.0,315.0,275.0,0.1,0.0,7.0,4.0,testThrottle,0.04,testTVCLower)

testRocket=Rocket(testPayload)
testRocket.addSerialStage(testUpperStage,0.5,2.0)
testRocket.addSerialStage(testLowerStage,0.0,2.0)

testPlace=np.array([0., 6371e3, 0.])
testVel=testEarth.absoluteWind(testPlace)

tea=np.linspace(-50,650,1000)
emm=np.zeros(len(tea))
thrust=np.zeros(len(tea))

for j in range(len(tea)):
    emm[j]=testRocket.mass(tea[j])
    thrust[j]=np.linalg.norm(testRocket.thrust(tea[j],testPlace,testEarth))
    
print("integrating...")
blah=testRocket.rk4Flight(testPlace,testVel,testEarth,0.0,800,True,7e6,1.5)
arblah=np.asarray(blah)
arlog=np.asarray(testRocket.flightlog)
print("done!")
print(testEarth.fiveKeplerianElements(arblah[-1,0:3],arblah[-1,3:6]))