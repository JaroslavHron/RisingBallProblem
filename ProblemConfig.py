from dolfin import *
from math import *
import os
import pickle
import mshr

class ProblemConfig:
    
    #####################################
    ##### Parameters of the problem #####
    # Specific numbers for the problem
    Re = 1.0
    T_spec = 1.0 # specific time

    # characteristic sizes
    x_c = 1.0
    t_c = 1.0
    v_c = 1.0

    # parameters of the ball
    alpha = 1.0
    R   = 1.0
    J_spec = 2.0/5.0  # specific moment of inertia, J = J_spec * m * R^2, J_sphere = 2.0/3.0, J_ball = 2.0/5.0
    m   = pi
    J   = 0.5 * pi
    m_e = pi
    J_e = 2.0/5.0 * pi

    # physical constants
    g   = 9.81  # normal g
    rho = 1e3  # density of water
    mu  = 1e-3  # dynamic viscosity of water
    nu  = 1e-6  # kinematic viscosity of water

    # Theoretical terminal velocities
    terminalStokesV = 1.0
    terminalOseenV = 1.0

    # Integration options

    dt = 0.1
    theta = 0.5
    n  = 10 # dt = T_spec/n
    N  = 100  # time of calculation will be N*T_spec
    T  = 10.0  # end time

    def __init__(self):
        self.Init(1e-3,0.9)
    
    ##
    # Direction of ball flow, 1=up, -1=down,
    # It is based on difference of the gravitational and buoyant force and the ball.
    def GetDir(self):
        if( self.alpha <= 1):
            return 'rising'
        else:
            return 'falling'

    def GetAcronym(self):
        return 'RisingBall_{0}_{1}'.format(self.R,self.x_c)

    def Save(self, fileName):
        file = open(fileName,'w')
        pickle.dump( self, file )        

    @staticmethod
    def Load(fileName):
        file = open(fileName,'r')
        return pickle.load( file )

    def Init(self, R, alpha, J_spec=2.0/5.0):

        # parameters of the ball
        self.alpha = alpha
        self.R = R
        self.J_spec = J_spec
        self.m = self.alpha * 4.0/3.0 * pi * self.R**3 * self.rho
        self.J = self.J_spec * self.m * self.R**2

        # Calculate terminal velocities
        self.terminalStokesV = 2.0/9.0 * abs(1-self.alpha) * self.R**2 * self.g / self.nu
        # self.terminalOseenV  =  ...

        # characteristic sizes
        self.x_c = self.R
        self.v_c = self.terminalStokesV
        self.t_c = self.x_c/self.v_c

        # effective mass and moment of inertia
        self.m_e = self.m * self.v_c / self.t_c
        self.J_e = self.J / self.t_c**2
        
        self.Re = self.v_c * self.x_c / self.nu
        self.T_spec = log(2) * self.m * self.v_c / (6.0 * pi * self.mu * self.x_c)
        self.dt = self.T_spec / self.n
        self.T  = self.T_spec * self.N

        self.CheckConsistency()
        
    def ChangePhysicalConstant(self, g=None, rho=None, mu=None):
        if g is not None:
            self.g   = g
        if rho is not None:
            self.rho = rho
        if mu is not None:
            self.mu  = mu
        self.nu  = self.mu/self.rho

        self.Init( self.R, self.alpha, self.J_spec )

    def ChangeTimeDiscretization(self, n=None, N=None):
        if n is not None:
            self.n = int(n)
        if N is not None:
            self.N = int(N)
        self.dt = self.T_spec / self.n
        self.T  = self.T_spec * self.N

        self.CheckConsistency()
        
    ##
    # Checks if this configuration is consistent i.e. checks that Re = x_c * v_c / nu etc.
    # @return True if it is consistent, raise exception AssertionError otherwise
    def CheckConsistency(self):

        def checkEqual( a, b ):
            if( abs((a-b)/a) < 1e-6 ):
                return True
            raise AssertionError('Configuration file is inconsistent, {0} != {1}'.format(a,b))

        checkEqual(self.m   , self.alpha * 4.0/3.0 * pi * self.R**3 * self.rho)
        checkEqual(self.J   , self.J_spec * self.m * self.R**2)

        checkEqual(self.x_c , self.R)
        checkEqual(self.v_c , self.x_c / self.t_c)
        
        checkEqual(self.m_e , self.m * self.v_c / self.t_c)
        checkEqual(self.J_e , self.J / self.t_c**2)

        checkEqual(self.Re  , self.v_c * self.x_c / self.nu)
        checkEqual(self.dt  , self.T_spec / self.n)
        checkEqual(self.T   , self.T_spec * self.N)

        return True
