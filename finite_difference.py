import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time as timel
### Initial conditions
nx=500 ## Number of grid points
nt=5000 ## Number of time steps
dx=1
dt=0.001
time2= np.linspace(0,nt*dt,nt+1)
p_old,p,p_new = (np.zeros(nx) for _ in range(3))
c = np.zeros(nx) + 700
# c= c + -300/249**2 * (np.arange(nx)-249)**2
domain =np.arange(0, nx*dx, dx)
### Source time function
f0=25.
t0=4./f0
time = np.linspace(0*dt, nt*dt, nt)

def src_func(time):
    return -2. * (time-t0) * (f0 ** 2) * (np.exp(-1.0 * (f0 **2) * (time - t0) ** 2))
isrc = round(nx/2)
### Boundary conditions
boundary = 'neumann'  

### Plot setup
fig,ax = plt.subplots(figsize=(10, 5))
line1,= ax.plot(domain, p, lw=2, color='blue')
# ax.set_xlim(4000,6000)
# five_pt_stencil(p, p_old, p_new, c, src, isrc, dt, dx, time2, nx, boundary)
# ### Solver
t=0
def five_pt_stencil(frame):
    global p,p_old,p_new, t 
    for x in np.arange(2, nx-2, 1):
        if x == isrc:
            p_new[x] = c[x]**2 * dt**2/ dx**2 *(-p[x-2]+16*p[x-1]-30*p[x]+16*p[x+1]-p[x+2])/12 + (src_func(t*dt)/dx * dt**2) +2*p[x] - p_old[x] 
        else:
            p_new[x] = c[x]**2 * dt**2/ dx**2 *(-p[x-2]+16*p[x-1]-30*p[x]+16*p[x+1]-p[x+2])/12 +2*p[x] - p_old[x] 
    p_old = p.copy()
    p=p_new.copy()
    if boundary == 'zero':
        p_new[0:1] = 0
        p_new[-2:] = 0
    elif boundary == 'neumann':
        p_new[0:1] = p_new[2]
        p_new[-2:] = p_new[-3]
    # elif boundary == 'absorbing':
    #     p_new[0:1] = p_new[2:]
    #     p_new[-2:] = 0
    line1.set_ydata(p)
    ax.set_ylim(np.min(p), np.max(p))
    t+=1
    if t%5==0:
        return(line1,)
ani = animation.FuncAnimation(fig, func= five_pt_stencil, frames=np.arange(0,500,5), interval=1, blit=False)
plt.show()