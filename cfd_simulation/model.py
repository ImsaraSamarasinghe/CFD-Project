import numpy as np
import matplotlib 
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
from matplotlib import cm

class Model:
    def __init__(self,lx,ly,nx,ny,nt):
        ## Dimensions
        self.lx = lx
        self.ly = ly
        self.nx = nx
        self.ny = ny
        self.dx = self.lx/(self.nx-1)
        self.dy = self.ly/(self.ny-1)
        self.nt = nt
        self.dt = 0.001

        self.x = None
        self.y = None
        self.X = None
        self.Y = None

        ## Physical constants
        self.rho = 1
        self.nu = 0.1
        self.c = 0.2
        self.dt = self.c*self.dx**2/self.nu
        ## Arrays
        self.u = np.zeros((self.ny, self.nx))
        self.v = np.zeros((self.ny, self.nx))
        self.p = np.zeros((self.ny, self.nx))
        self.b = np.zeros((self.ny, self.nx))
    
    def CreateGrid(self):
        self.x = np.linspace(0.0,self.lx,self.nx)
        self.y = np.linspace(0.0,self.ly,self.ny)
        self.X, self.Y = np.meshgrid(self.x,self.y)

        # plot the mesh
        plt.plot(self.X,self.Y,marker='o',linestyle='',color='k')
        plt.savefig('mesh.png')
    
    def build_up_b(self):
    
        self.b[1:-1, 1:-1] = (self.rho * (1 / self.dt * 
                        ((self.u[1:-1, 2:] - self.u[1:-1, 0:-2]) / 
                        (2 * self.dx) + (self.v[2:, 1:-1] - self.v[0:-2, 1:-1]) / (2 * self.dy)) -
                        ((self.u[1:-1, 2:] - self.u[1:-1, 0:-2]) / (2 * self.dx))**2 -
                        2 * ((self.u[2:, 1:-1] - self.u[0:-2, 1:-1]) / (2 * self.dy) *
                            (self.v[1:-1, 2:] - self.v[1:-1, 0:-2]) / (2 * self.dx))-
                            ((self.v[2:, 1:-1] - self.v[0:-2, 1:-1]) / (2 * self.dy))**2))
    
    def pressure_poisson(self):
        pn = np.empty_like(self.p)
        pn = self.p.copy()
        
        for q in range(self.nt//2):
            pn = self.p.copy()
            self.p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * self.dy**2 + 
                            (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * self.dx**2) /
                            (2 * (self.dx**2 + self.dy**2)) -
                            self.dx**2 * self.dy**2 / (2 * (self.dx**2 + self.dy**2)) * 
                            self.b[1:-1,1:-1])

            self.p[:, -1] = self.p[:, -2] # dp/dx = 0 at x = 2
            self.p[0, :] = self.p[1, :]   # dp/dy = 0 at y = 0
            self.p[:, 0] = self.p[:, 1]   # dp/dx = 0 at x = 0
            self.p[-1, :] = 0        # p = 0 at y = 2
    
    def cavity_flow(self):
        un = np.empty_like(self.u)
        vn = np.empty_like(self.v)
        
        for n in range(self.nt):
            un = self.u.copy()
            vn = self.v.copy()
            
            self.build_up_b()
            self.pressure_poisson()
            
            self.u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
                            un[1:-1, 1:-1] * self.dt / self.dx *
                            (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                            vn[1:-1, 1:-1] * self.dt / self.dy *
                            (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                            self.dt / (2 * self.rho * self.dx) * (self.p[1:-1, 2:] - self.p[1:-1, 0:-2]) +
                            self.nu * (self.dt / self.dx**2 *
                            (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                            self.dt / self.dy**2 *
                            (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

            self.v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
                            un[1:-1, 1:-1] * self.dt / self.dx *
                        (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                            vn[1:-1, 1:-1] * self.dt / self.dy *
                        (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                            self.dt / (2 * self.rho * self.dy) * (self.p[2:, 1:-1] - self.p[0:-2, 1:-1]) +
                            self.nu * (self.dt / self.dx**2 *
                        (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                            self.dt / self.dy**2 *
                        (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

            self.u[0, :]  = 0
            self.u[:, 0]  = 0
            self.u[:, -1] = 0
            self.u[-1, :] = 1    # set velocity on cavity lid equal to 1
            self.v[0, :]  = 0
            self.v[-1, :] = 0
            self.v[:, 0]  = 0
            self.v[:, -1] = 0

    def plot(self):

        fig = plt.figure(figsize=(11, 7), dpi=100)
        plt.contourf(self.X, self.Y, self.p, alpha=0.5, cmap=cm.viridis)
        plt.colorbar()
        plt.contour(self.X, self.Y, self.p, cmap=cm.viridis)
        plt.streamplot(self.X, self.Y, self.u, self.v)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.savefig('FlowField.png')