from dolfin import *
import os

class RisingBallProblem:

    # Configuration file
    mainCfg = []
    meshCfg = []

    # Mesh
    mesh    = []
    bndry   = []
    dS      = []
    
    # Function spaces
    V = []#VectorFunctionSpace(mesh, "CG", 2)
    P = []#FunctionSpace(mesh, "CG", 1)
    R = []#FunctionSpace(mesh, "R", 0)
    WD = []#MixedFunctionSpace([V, P, R, R, R, R, R, R])  # Space for dynamic problem
    WS = []#MixedFunctionSpace([V, P])                    # Space for static problem

    # Solution of static problem
    ws = []
    # Solution of dynamic problem
    wd = []
    wd0= []    

    # Boundary condition
    # BC for dynamic problem
    # Boundary condition base functions
    QX = []#Function(V)
    QY = []#Function(V)
    QZ = []#Function(V)
    WX = []#Function(V)
    WY = []#Function(V)
    WZ = []#Function(V)
    # BC for static problem
    # Expression used to define boundary condition on the ball in the static problem
    bc_vb = Expression(("qx+wy*x[2]-wz*x[1]","qy+wz*x[0]-wx*x[2]","qz+wx*x[1]-wy*x[0]"),qx=0.0,qy=0.0,qz=0.0,wx=0.0,wy=0.0,wz=0.0)
    # This is speed of the ball used in convection term in the static problem
    vb    = Expression(("qx","qy","qz"),qx=0.0,qy=0.0,qz=0.0)

    timeStep = Expression("dt",dt=1.0)

    # Solvers
    # Solvers used to create solution in next time step
    NavierStokesSolver = []
    NavierStokesProblem = []
    NavierStokesStaticSolver = []
    NavierStokesStaticProblem = []
    StokesStaticSolver = []
    StokesStaticProblem = []
    OseenStaticSolver = []
    OseenStaticProblem = []

    ###################
    ### Constructor ###
    def __init__(self,mainConfig,meshConfig):
        self.mainCfg = mainConfig
        self.meshCfg = meshConfig

        #parameters['form_compiler']['representation'] = 'uflacs'
        parameters['form_compiler']['optimize'] = True
        parameters['form_compiler']['quadrature_degree'] = 4

        self.mesh,self.bndry = self.meshCfg.LoadH5Mesh()

        self.dS = Measure("ds",domain=self.mesh,subdomain_data=self.bndry)
        self.timeStep.dt = self.mainCfg.dt

        self.V = VectorFunctionSpace(self.mesh, "CG", 2)
        self.P = FunctionSpace(self.mesh, "CG", 1)
        self.R = FunctionSpace(self.mesh, "R", 0)
        self.WD = MixedFunctionSpace([self.V, self.P, self.R, self.R, self.R, self.R, self.R, self.R])  # Space for dynamic problem
        self.WS = MixedFunctionSpace([self.V, self.P])                    # Space for static problem
        self.ws = Function(self.WS)
        self.wd = Function(self.WD)
        self.wd0= Function(self.WD)
        self.QX = Function(self.V)
        self.QY = Function(self.V)
        self.QZ = Function(self.V)
        self.WX = Function(self.V)
        self.WY = Function(self.V)
        self.WZ = Function(self.V)
        
        self.InitStaticSolvers()
        self.InitBCBase()
        self.InitDynamicSolvers()

    #######################################
    ### Solver initialization functions ###
    def InitStaticSolvers(self):
        v_,p_ = TestFunctions(self.WS)
        vt,pt = TrialFunctions(self.WS)

        bcs = self.InitStaticBC()

        v,p = split(self.ws)
        SF  = self.SA(vt,pt,v_) + self.PA(vt,p_)
        NSF = self.SNSA(v,p ,v_) + self.PA(v ,p_)
        self.StokesStaticSolver,self.StokesStaticProblem = self.InitLinearSolver( SF, self.ws, bcs )
        self.NavierStokesStaticSolver,self.NavierStokesStaticProblem = self.InitNonlinearSolver( NSF, self.ws, bcs )
        
    def InitDynamicSolvers(self):
        v_,p_,qx_,qy_,qz_,wx_,wy_,wz_ = TestFunctions(self.WD)
        v ,p ,qx ,qy ,qz ,wx ,wy ,wz  = split(self.wd )
        v0,p0,qx0,qy0,qz0,wx0,wy0,wz0 = split(self.wd0)

        bcs = self.InitDynamicBC()

        cfg = self.mainCfg
        domainNorm = Constant( 1.0/(4*pi) ) # The ODE's are integrated over the surface of the ball, so we need to renormalize them by ball's surface area
        u  = self.GetFullVelocity(self.wd )
        u0 = self.GetFullVelocity(self.wd0)
        #
        
        F =  inner(u-u0,v_)/self.timeStep * dx + domainNorm * ( (qx-qx0)*qx_ + (qy-qy0)*qy_ + (qz-qz0)*qz_ + (wx-wx0)*wx_ + (wy-wy0)*wy_ + (wz-wz0)*wz_)/self.timeStep * self.dS(3) + \
             (1-cfg.theta) * (self.DNSA(u0,p,qx0,qy0,qz0,v_) + self.PA(v,p_)  - self.BFA(u0,p ,qx_,qy_,qz_) - self.BMA(u0,p ,wx_,wy_,wz_)) + \
             (  cfg.theta) * (self.DNSA(u ,p,qx ,qy ,qz ,v_) + self.PA(v,p_)  - self.BFA(u ,p ,qx_,qy_,qz_) - self.BMA(u ,p ,wx_,wy_,wz_))
        # note : BC base is divergence free, thus v in pressure term is sufficient

        self.NavierStokesSolver,self.NavierStokesProblem = self.InitNonlinearSolver( F, self.wd, bcs )
        

    def InitLinearSolver(self, F, w, bcs ):
        problem = LinearVariationalProblem(lhs(F), rhs(F), w, bcs)
        solver  = LinearVariationalSolver(problem)
        
        prm = solver.parameters
        prm['linear_solver'] = 'mumps'
        if False:
            prm["linear_solver"] = "gmres"
            prm["preconditioner"] = "amg"
            prm["krylov_solver"]["absolute_tolerance"] = 1E-6
            prm["krylov_solver"]["relative_tolerance"] = 1E-5
            prm["krylov_solver"]["maximum_iterations"] = 1000
            prm["krylov_solver"]["gmres"]["restart"] = 40
            prm["krylov_solver"]["monitor_convergence"] = True

        return solver,problem
        
    def InitNonlinearSolver(self, F, w, bcs ):
        J = derivative( F, w )
        problem  = NonlinearVariationalProblem( F, w, bcs, J )
        solver = NonlinearVariationalSolver( problem )

        prm = solver.parameters
        #info(prm,True)
        prm['nonlinear_solver'] = 'newton'
        prm['newton_solver']['absolute_tolerance'] = 1E-6
        prm['newton_solver']['relative_tolerance'] = 1e-5
        prm['newton_solver']['maximum_iterations'] = 25
        prm['newton_solver']['relaxation_parameter'] = 1.0
        prm['newton_solver']['linear_solver'] = 'mumps'
        prm['newton_solver']['report'] = True
        if False:
            prm['newton_solver']["linear_solver"] = "gmres"
            prm['newton_solver']["preconditioner"] = "bjacobi"
            prm['newton_solver']["krylov_solver"]["absolute_tolerance"] = 1E-9
            prm['newton_solver']["krylov_solver"]["relative_tolerance"] = 1E-7
            prm['newton_solver']["krylov_solver"]["maximum_iterations"] = 300
            prm['newton_solver']["krylov_solver"]["gmres"]["restart"] = 40
            prm['newton_solver']["krylov_solver"]["nonzero_initial_guess"] = True
            prm['newton_solver']["krylov_solver"]["monitor_convergence"] = True
    
        return solver,problem

    def InitStaticBC(self):
        noslip = Constant((0.0, 0.0, 0.0))
        bcv_walls   = DirichletBC(self.WS.sub(0), noslip, self.bndry, 1)
        #bcv_outflow = DirichletBC(self.WS.sub(0), noslip, self.bndry, 2)
        bcv_ball    = DirichletBC(self.WS.sub(0), self.bc_vb, self.bndry, 3)

        # Collect boundary conditions
        bcs = [bcv_ball,bcv_walls]
        return bcs
        
    def InitDynamicBC(self):
        noslip = Constant((0.0, 0.0, 0.0))
        bcv_walls   = DirichletBC(self.WS.sub(0), noslip, self.bndry, 1)
        #bcv_outflow = DirichletBC(self.WS.sub(0), noslip, self.bndry, 2)
        bcv_ball    = DirichletBC(self.WS.sub(0), noslip, self.bndry, 3)

        # Collect boundary conditions
        bcs = [bcv_ball,bcv_walls]
        return bcs

    def InitBCBase(self,force=False):

        bfNames = ['QX','QY','QZ','WX','WY','WZ']
        bfFilePath   = os.path.join(self.meshCfg.meshDirectory,'bcBaseFunctions',self.meshCfg.GetAcronym()+'_bcBaseFunctions.h5')
        xdmfFilePath = os.path.join(self.meshCfg.meshDirectory,'bcBaseFunctions',self.meshCfg.GetAcronym()+'_bcBaseFunctions.xdmf')

        bfLoaded = False
        # Load base functions if they were calculated before
        if os.path.isfile(bfFilePath) and force==False:
            print('Loading bc base functions.')
            bfFile = HDF5File(mpi_comm_world(),bfFilePath,'r')
            
            i = 0
            try:
                for BF in [self.QX,self.QY,self.QZ,self.WX,self.WY,self.WZ]:
                    bfFile.read(BF,bfNames[i])
                    i += 1
                bfLoaded = True
                print('Loading done.')
            except:
                print('Could not load bc base functions!')
                bfLoaded = False


        # Calculate base functions
        if not bfLoaded:
            print('Calculating bc base functions.')
            if os.path.isfile(bfFilePath): # HDF5File has a problem with rewriting existing files for some reason, so we delete manualy
                os.remove(bfFilePath)

            bfFile   = HDF5File(mpi_comm_world(),bfFilePath,'w')
            #xdmfFile = XDMFFile(mpi_comm_world(),xdmfFilePath)
            
            i = 0
            for BF in [self.QX,self.QY,self.QZ,self.WX,self.WY,self.WZ]:
                
                print('Calculating base function: {0}'.format(bfNames[i]))
                fa = FunctionAssigner( self.V, self.WS.sub(0) )
                Q = [0,0,0,0,0,0]
                Q[i] = 1.0
                
                self.SetStaticBC(Q)
                self.StokesStaticSolver.solve()
                fa.assign(BF,self.ws.sub(0))
                bfFile.write(BF,bfNames[i])

                # It does not work for some reason
                # Save base functions for viewing
                #v,p = self.ws.split()
                #v.rename('v','velocity_'+bfNames[i])
                #p.rename('P'+bfNames[i],'pressure_'+bfNames[i])
                #xdmfFile << v
                #xdmfFile << p
                
                i += 1
            #self.ws.sub(0).rename('v','velocity')
            #self.ws.sub(1).rename('p','pressure')
        
    def SetStaticBC(self,Q):
        self.bc_vb.qx = Q[0]
        self.bc_vb.qy = Q[1]
        self.bc_vb.qz = Q[2]
        self.bc_vb.wx = Q[3]
        self.bc_vb.wy = Q[4]
        self.bc_vb.wz = Q[5]
        self.vb.qx = Q[0]
        self.vb.qy = Q[1]
        self.vb.qz = Q[2]

    def GetFullVelocity(self, w):
        v,p,qx,qy,qz,wx,wy,wz = split(w)
        return v + qx*self.QX + qy*self.QY + qz*self.QZ + wx*self.WX + wy*self.WY + wz*self.WZ
        
    ##############################
    ### Weak formulation forms ###

    # Viscosity form
    def DA(self,v,_v):
        cfg = self.mainCfg
        D   = sym(grad(v))
        Re  = Constant( cfg.Re )
        return 2.0/Re * inner(D, grad(_v)) * dx
        #return inner(D, grad(_v)) * dx

    # Pressure from
    # @param v velocity
    def PA(self,v,p):
        cfg = self.mainCfg
        return inner(div(v),p) * dx

    # Convection form
    # @param v velocity to be convected
    # @param u convection velocity
    def CA(self,v,u,_v):
        return inner(grad(v)*u,_v) * dx

    # Stokes form
    def SA(self,v,p,_v):
        return - self.PA(_v,p) + self.DA(v,_v)

    # Oseen form
    def OA(self,v,p,_v):
        return self.SA(v,p,_v) - self.CA(v,self.vb,_v) 

    # Navier Stokes form static
    def SNSA(self,v,p,_v):
        return self.OA(v,p,_v) + self.CA(v,v,_v)

    def DNSA(self,v,p,qx,qy,qz,_v):
        U = qx*Constant((1,0,0)) + qy*Constant((0,1,0)) + qz*Constant((0,0,1))
        return self.SA(v,p,_v) - self.CA(v,U,_v) + self.CA(v,v,_v)

    # Dynamic ODE of ball
    # @param u full fluid sped, i.e. u = v + qx*QX + qy*QY + ...
    def BFA(self,v,p,_qx,_qy,_qz):
        ns = Constant(1.0/(4.0*pi)) 
        F  = Constant(ns*(self.F_g()+self.F_vz()))
        dS = Measure("ds",subdomain_data=self.bndry)
        return (_qx*self.F_f(v,p,0) + _qy*self.F_f(v,p,1) + _qz*(self.F_f(v,p,2)+F)) * self.dS(3)

    def BMA(self,v,p,_wx,_wy,_wz):
        return (_wx*self.M_f(v,p,0) + _wy*self.M_f(v,p,1) + _wz*self.M_f(v,p,2) ) * self.dS(3)


    #######################
    ### Force functions ###
    def F_g(self):
        cfg = self.mainCfg
        return -cfg.x_c * cfg.g / cfg.v_c**2

    def F_vz(self):
        cfg = self.mainCfg
        return cfg.x_c * cfg.g / (cfg.alpha * cfg.v_c**2)
    
    def F_f(self,v,p,i):
        cfg = self.mainCfg
        n = -FacetNormal(self.mesh) # outer normal to mesh is inner normal of ball
        D = sym(grad(v))
        Dn = D*n
        #C1 = Constant(      cfg.x_c**2 * cfg.v_c**2 * cfg.rho / cfg.m_e)
        #C2 = Constant(2.0 * cfg.mu     * cfg.v_c    * cfg.x_c / cfg.m_e)
        C1 = Constant(3.0 / (4.0 * pi * cfg.alpha ))
        C2 = Constant(6.0 / (4.0 * pi * cfg.alpha * cfg.Re))
        F = (- C1 * p * n[i]  + C2 * Dn[i] )
        return F

    def M_f(self,v,p,i):
        cfg = self.mainCfg
        n = -FacetNormal(self.mesh) # outer normal to mesh is inner normal of ball
        D = sym(grad(v))
        Dn = D*n
        x = Expression(("x[0]","x[1]","x[2]"))
        print('Moment se nepocita spravne!! musim byt vyjadren pomoci J_spec')
        C1 = Constant( 6.0 / (4.0 * pi * cfg.alpha)) 
        C2 = Constant(12.0 / (4.0 * pi * cfg.alpha))
        j = ( i + 1 ) % 3
        k = ( i + 2 ) % 3
        M = (- C1 * p * (x[j]*n[k]-x[k]*n[j]) + C2 * (x[j]*Dn[k]-x[k]*Dn[j]) )
        return M

    # Total forces and moments on the ball
    def ForceOnBall(self,v,p,i):
        F = self.F_f(v,p,i)
        if i==2:
            ns = Constant(1.0/(4.0*pi)) # Integration normalization
            F += Constant(ns*(self.F_g()+self.F_vz()))
        return F

    def TotalForceOnBall(self):
        cfg = self.mainCfg
        u = self.ws.sub(0)
        p = self.ws.sub(1)
        Fx = assemble( self.ForceOnBall(u, p, 0) * self.dS(3) ) * cfg.m_e
        Fy = assemble( self.ForceOnBall(u, p, 1) * self.dS(3) ) * cfg.m_e
        Fz = assemble( self.ForceOnBall(u, p, 2) * self.dS(3) ) * cfg.m_e
        return [Fx,Fy,Fz]

    # Only force from fluid
    # @return Force [Fx,Fy,Fz] in units 
    def FluidForceOnBall(self,v,p):
        cfg = self.mainCfg
        Fx = assemble( self.F_f(v,p,0)*self.dS(3) ) * cfg.m_e
        Fy = assemble( self.F_f(v,p,1)*self.dS(3) ) * cfg.m_e
        Fz = assemble( self.F_f(v,p,2)*self.dS(3) ) * cfg.m_e
        return [Fx,Fy,Fz]

    # Moment force on ball from fluid
    # @return Moment force [Mx,My,Mz] in units
    def FluidMomentOnBall(self,v,p):
        cfg = self.mainCfg
        Mx = assemble( self.M_f(v,p,0)*self.dS(3) ) * cfg.J_e
        My = assemble( self.M_f(v,p,1)*self.dS(3) ) * cfg.J_e
        Mz = assemble( self.M_f(v,p,2)*self.dS(3) ) * cfg.J_e
        return [Mx,My,Mz]
