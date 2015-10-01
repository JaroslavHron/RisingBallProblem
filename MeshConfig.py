import os
import os.path
import dolfin
from math import *

##
# Class containing mesh parameters and has functions to generate such a mesh in 2d or 3d
#
# Mesh Shape:
#    Mesh is always box with cut out unit ball in the origin.
#    Most of the parameters affect the density of the mesh.
#    
class MeshConfig:

    #######################
    ### Mesh Parameters ###
    # Box size
    L1 = 30.0
    L2 = 20.0
    # Element size
    s = 0.4
    a = 0.4 #0.4
    b = 0.05
    alpha = 1.0
    # Shape parameters
    d1    = 2.0
    d2    = 10.0
    theta = pi/2.0
    # Dimension
    dim = 3

    # Directory to store meshes
    meshDirectory = ""
    templateGeoFile = ""

    #####################
    #### Constructor ####

    ##
    # Construct MeshConfig with default parameters
    # @param directory Directory for meshes, it must contain template geo file "boxWithHole1.geo".
    #        If the mesh directory is in current working directory, than than the argument can be left as default
    def __init__(self,directory=None):
        if(directory is None):
            directory = os.path.join(os.getcwd(),'mesh')        

        self.templateGeoFile = os.path.join(directory,'template','boxWithHole1.geo')

        if(not os.path.isfile(self.templateGeoFile)):
            raise ValueError('File {0} could not be found!'.format(self.templateGeoFile))
        self.meshDirectory = directory
    
    #########################
    ### Utility functions ###
    # Acronym of this mesh
    def GetAcronym(self):
        return 'box_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}_{9}'.format( int(self.dim), float(self.L1), float(self.L2), float(self.s), float(self.a), float(self.b), float(self.alpha), float(self.d1), float(self.d2), float(self.theta))

    #def LoadFromAcronym( word ):
    # Needs implementation
        
    ##
    # Get mesh file path
    # @param suffix Can be '.geo' or '.msh' or '.xml'
    def GetFilePath(self,suffix):
        return os.path.join(self.meshDirectory,self.GetAcronym()+suffix)
    
    ##############################
    ### Construction functions ###
    ##
    # Creates geo file with this configuration
    def GenerateGeo(self):

        filePath = self.GetFilePath('.geo')
        geoFile = open(filePath,'w')

        if(self.a<self.b):
            raise ValueError("Invalid mesh settings: a can not be smaller than b")

        geoFile.write('L1 = {0};\n'.format(self.L1))
        geoFile.write('L2 = {0};\n'.format(self.L2))
        geoFile.write('s  = {0};\n'.format(self.s))
        geoFile.write('a  = {0};\n'.format(self.a-self.b))
        geoFile.write('b  = {0};\n'.format(self.b))
        geoFile.write('alpha = {0};\n'.format(self.alpha))
        geoFile.write('d1  = {0};\n'.format(self.d1))
        geoFile.write('d2  = {0};\n'.format(self.d2))
        geoFile.write('theta = {0};\n'.format(self.theta))
        
        geoFile.close()

        os.system('cat {0} >> {1}'.format(self.templateGeoFile,filePath))

    def GenerateMsh(self,force=False):
        geoFilePath = self.GetFilePath('.geo')
        mshFilePath = self.GetFilePath('.msh')
        if(os.path.isfile(geoFilePath) and \
           (not os.path.isfile(mshFilePath) or force==True)):
            os.system('gmsh -3 -algo del3d -o {0} {1}'.format(mshFilePath,geoFilePath))
        else:
            if(not os.path.isfile(geoFilePath)):
                raise IOError('Could not find geo file {0}!'.format(geoFilePath))

    def GenerateXml(self,force=False):
        mshFilePath = self.GetFilePath('.msh')
        xmlFilePath = self.GetFilePath('.xml')
        if(os.path.isfile(mshFilePath) and  \
           (not os.path.isfile(xmlFilePath) or force==True)):
            os.system('dolfin-convert {0} {1}'.format(mshFilePath,xmlFilePath))
        else:
            if(not os.path.isfile(mshFilePath)):
                raise IOError('Could not find msh file {0}!'.format(mshFilePath))

    def GenerateH5(self,force=False):
        h5FilePath  = self.GetFilePath('.h5')

        if(not os.path.isfile(h5FilePath) or \
           force==True):
            h5File = dolfin.HDF5File(dolfin.mpi_comm_world(),h5FilePath,'w')
            mesh,bndry = self.LoadXMLMesh(force)        
            h5File.write(mesh,'mesh')
            h5File.write(bndry,'bndry')
                
    def GenerateMesh(self,force=False):
        self.GenerateGeo()
        self.GenerateMsh(force)
        self.GenerateXml(force)
            
    def LoadXMLMesh(self,force=False):
        # Make sure that mesh exists
        self.GenerateMesh(force)
        
        xmlFilePath = self.GetFilePath('.xml')
        bndryFilePath = self.GetFilePath('_facet_region.xml')
        mesh = []
        bndry = []
        if(os.path.isfile(xmlFilePath)):
            mesh = dolfin.Mesh(xmlFilePath)
        else:
            raise IOError('Could not find xml file {0}!'.format(xmlFilePath))
            
        if(os.path.isfile(bndryFilePath)):
            bndry = dolfin.MeshFunction("size_t",mesh, bndryFilePath)
        else:
            raise IOError('Cound not find xml file {0}!'.format(bndryFilePath))
        
        return [mesh,bndry]
        
    def LoadH5Mesh(self,force=False):
        # Make sure that mesh exists
        self.GenerateH5(force)

        mesh  = dolfin.Mesh()
        bndry = []

        h5FilePath  = self.GetFilePath('.h5')
        if(os.path.isfile(h5FilePath)):
            h5File = dolfin.HDF5File(dolfin.mpi_comm_world(),h5FilePath,'r')
            h5File.read(mesh,'mesh',False)
            a = dolfin.FacetFunction("size_t", mesh)
            h5File.read(a,'bndry')
            bndry = a
        else:
            raise IOError('Could not find h5 file {0}!'.format(h5FilePath))
            
        return [mesh,bndry]


    def ViewInGmsh(self):
        mshFilePath = self.GetFilePath('.msh')
        posFilePath = os.path.join(self.meshDirectory,'gmshCfg.pos')
        if(not os.path.isfile(posFilePath)):
            posFilePath = ""
        if(os.path.isfile(mshFilePath)):
            os.system('gmsh {0} {1}'.format(mshFilePath,posFilePath))
        else:
            raise IOError('Could not find msh file {0}!'.format(mshFilePath))

