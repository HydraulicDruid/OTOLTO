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

eccentricityWeightCap=10000
inclinationWeightCap=1000
excessFuelWeightingFactor=200 #per kilogram - 200 is unusually high, of course

desiredEccentricity=0.0
desiredInclination=0.0

testAtmosphere=Atmosphere(0,1.2,29.26,260)
testEarth=Planet(5.97219e24, 6371e3, 86164.1, testAtmosphere)
coords=np.array([0.0, 6371010.0, 0.0])
velocity=np.array([0.0,0.0,0.0])

#print(testEarth.absoluteWind(coords))
#print(str(testEarth.get_atm_rho(coords)))
#print(testEarth.fiveKeplerianElements(coords,velocity))

testUpperStage=Stage(120.0,7.0,350.0,0.0,0.7,0.7,0.28,0.3,0.3)
testLowerStage=Stage(1400.0,7.2,290.0,265.0,0.7,0.8,0.28,7.0,4.0)

#testSSTO=Stage(1.25e5,25,450.0,350.0,0.8,0.8,12.6,400.0,20.0)

propellantSchedule=np.array([[0.0,175.0,360.0,3e7],[7.0,5.0,0.3,0.3]])

#rows are time, theta, phi (theta is azimuth and phi is polar)
tvcSchedule=np.array([[0.0,10.1,48.5,189.8,365.6,3e7], [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2], [pi/2+0.0,pi/2+0.0,pi/2+0.204,pi/2+0.870,pi/2+1.468,pi/2+1.468]])

testRocket=Rocket(propellantSchedule,tvcSchedule)
testRocket.addStage(testUpperStage,1,5)
testRocket.addStage(testLowerStage,20,0)

#testRocket.addStage(testSSTO, 3e7, 0)

blah=testRocket.simulateFlight(coords,velocity,testEarth,600.0,8.0,True,7e6)
arblah=np.asarray(blah).T

esT=np.linspace(-pi,pi,2000)
esX,esY=6371e3*np.sin(esT),6371e3*np.cos(esT)
pylab.plot(esX,esY)

xwidth=max(arblah[1])-min(arblah[1])
yheight=max(arblah[2])-min(arblah[2])

pylab.xlim( min(arblah[1])-0.1*xwidth, max(arblah[1])+0.1*xwidth )
pylab.ylim( min(arblah[2])-0.1*yheight, max(arblah[2])+0.1*yheight )
pylab.plot(arblah[1],arblah[2])

kepEls=testEarth.fiveKeplerianElements(arblah[1:4,-1],arblah[4:7,-1])

print(kepEls)
print(arblah[17,-1])

fitness=arblah[17, -1]*excessFuelWeightingFactor + min(eccentricityWeightCap, 1/abs(kepEls[0]-desiredEccentricity)) + min(inclinationWeightCap, 1/abs(kepEls[2]-desiredInclination))
print(fitness)