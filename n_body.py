import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib.animation import FuncAnimation

class Particle:
    def __init__(self,r_0,v_0,m):
        self.r_0=r_0
        self.v_0=v_0
        self.r=r_0
        self.v=v_0
        self.m=m
    def update(self,F,dt):
        a=F/self.m
        self.v += a*dt
        self.r += self.v*dt

class Simulation:
    def __init__(self,boxsize,dt,particles=None,particle_num=None):
        self.boxsize=boxsize
        self.dt=dt
        self.particles=particles
        self.particle_num=len(particles)
        if particles and particle_num:
            raise ValueError("Only one of particles or particle_num should be defined")
        if not particles and not particle_num:
            raise ValueError("At least one of particles or particle_num should be defined")
    def get_forcematrix(self):
        Fij= np.zeros((self.particle_num,self.particle_num,3))
        for j in range(1,self.particle_num):
            for i in range(j):
                Fij[i,j,:]=self.get_gravity(self.particles[i],self.particles[j])
                Fij[j,i,:]=-Fij[i,j,:]
        return Fij
    def get_gravity(self,p1,p2):
        dr=p1.r-p2.r
        mag_dr=np.sqrt(np.dot(dr,dr))
        return -p1.m*p2.m/mag_dr**3*dr
    def get_forces(self):
        farray=np.zeros((self.particle_num,3))
        fij=self.get_forcematrix()
        for i in range(self.particle_num):
            farray=fij.sum(axis=1)
        return farray
    def step(self):
        farray=self.get_forces()
        for i in range(self.particle_num):
            self.particles[i].update(farray[i],self.dt)
        
M=10000
m=1
r=10
p1=Particle(m=M,r_0=np.zeros(3),v_0=np.zeros(3))
p2=Particle(m=m,r_0=np.array([r,0.0,0.0]),v_0=np.array([0.0,np.sqrt(M/r),0.0]))

simulation=Simulation(1,1e-4,[p1,p2])

nsteps=63000

pos2=np.zeros((nsteps,3))
for i in range(nsteps):
    pos2[i]=p2.r
    simulation.step()

Fig, ax =plt.subplots()

x1data,y1data=[],[]
x2data,y2data=[],[]

plot1, =plt.plot([],[],"bo")
plot2, =plt.plot([],[],"ro")

def init():
    ax.set_xlim(-11,11)
    ax.set_ylim(-11,11)
    ax.set_aspect("equal")

    return plot1,plot2

def update(frame):
    # x1data.append(pos2[frame,0])
    # y1data.append(pos2[frame,1])
    x1data = pos2[frame,0]
    y1data = pos2[frame,1]

    x2data.append(0)
    y2data.append(0)

    plot1.set_data(x1data,y1data)
    plot2.set_data(x2data,y2data)
    return plot1,plot2

animation=FuncAnimation(Fig,update,frames=[i*200 for i in range(nsteps//200)],blit=True, init_func=init)

animation.save("1planet.mp4",fps=60)
plt.close()
M=100000
m=1
m_moon=0.01
r=10
r_moon=10.01
p1=Particle(m=M,r_0=np.zeros(3),v_0=np.zeros(3))
p2=Particle(m=m,r_0=np.array([r,0.0,0.0]),v_0=np.array([0.0,np.sqrt(M/r),0.0]))
p3=Particle(m=m_moon,r_0=np.array([r_moon,0.0,0.0]),v_0=np.array([0.0,np.sqrt(m/(r_moon-r))+np.sqrt(M/r_moon),0.0]))
simulation=Simulation(1,1e-5,[p1,p2,p3])
nsteps=62831
pos2=np.zeros((nsteps,3))
pos3=np.zeros((nsteps,3))
for i in range(nsteps):
    pos2[i]=p2.r
    pos3[i]=p3.r
    simulation.step()
Fig, ax =plt.subplots()

x1data,y1data=[],[]
x2data,y2data=[],[]
x3data,y3data=[],[]

plot1, =plt.plot([],[],"bo",label="Earth",mec="b")
plot2, =plt.plot([],[],"ro")
plot3, =plt.plot([],[],"go", label="Moon",markersize=5, mec="g")

def init():
    ax.set_xlim(-1.3,1.3)
    ax.set_ylim(9.8,10.1)
    plt.legend()
    # ax.set_aspect("equal")

    return plot1,plot2,plot3

def update(frame):
    # x1data.append(pos2[frame,0])
    # y1data.append(pos2[frame,1])
    x1data = pos2[frame,0]
    y1data = pos2[frame,1]

    x2data.append(0)
    y2data.append(0)

    x3data = pos3[frame,0]
    y3data = pos3[frame,1]

    plot1.set_data(x1data,y1data)
    plot2.set_data(x2data,y2data)
    plot3.set_data(x3data,y3data)

    return plot1,plot2,plot3

animation=FuncAnimation(Fig,update,frames=[i for i in range(14_300,17_000, 5)],blit=True, init_func=init)

animation.save("earth_moon.mp4",fps=90)
plt.close()