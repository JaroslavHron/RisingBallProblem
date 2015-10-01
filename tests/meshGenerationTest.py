import os
import sys
sys.path.append('../')
from MeshConfig import *
from dolfin import *

cfg = MeshConfig('../mesh/')

cfg.L1 = 40
cfg.L2 = 30
cfg.d2 = 20.0
cfg.theta = 0.8*pi/2
cfg.s = 0.4

cfg.GenerateGeo()
cfg.GenerateMsh(False)
cfg.GenerateMesh(False)
cfg.ViewInGmsh()

#mesh,bndry = cfg.LoadH5Mesh()
#plot(bndry,wireframe=True,interactive=True)

