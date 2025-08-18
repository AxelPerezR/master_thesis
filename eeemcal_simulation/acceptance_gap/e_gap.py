from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV, deg
SIM = DD4hepSimulation() 

SIM.gun.thetaMin = 90*deg
SIM.gun.thetaMax = 178*deg
SIM.gun.distribution = "eta"
SIM.gun.momentumMin = 20*GeV
SIM.gun.momentumMax = 20*GeV
SIM.gun.particle = "e-"
SIM.enableGun = True
SIM.numberOfEvents = 10000
SIM.outputFile = "e_20GeV_Gap_200cm_90to178deg_1e4.edm4hep.root"