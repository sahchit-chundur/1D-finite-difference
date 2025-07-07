import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import writers
### Initial conditions
nx=500 ## Number of grid points
nt=5000 ## Number of time steps
dx=1
dt=0.001
time2= np.linspace(0,nt*dt,nt+1)
p_old,p,p_new = (np.zeros(nx) for _ in range(3))
c0=700
c = np.zeros(nx) + c0
# c[0:int(nx/4)] = c0 + 200
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
boundary = 'zero'  

### Plot setup

cfl = c0*dt/dx
if cfl > 1:
    print(f'Warning: CFL condition violated! CFL = {cfl:.2f} > 1.0')
else:
    print(f'CFL condition satisfied. CFL = {cfl:.2f} <= 1.0')
fig,ax = plt.subplots(figsize=(10, 5))
line1,= ax.plot(domain, p, lw=2, color='blue')
# ### Solver
t=0
def five_pt_stencil(frame):
    global p,p_old,p_new, t 
    for x in np.arange(2, nx-2, 1):
        if x == isrc:
            p_new[x] = c[x]**2 * dt**2/ dx**2 *(-p[x-2]+16*p[x-1]-30*p[x]+16*p[x+1]-p[x+2])/12 + (src_func(t*dt)/dx * dt**2) +2*p[x] - p_old[x] 
        else:
            p_new[x] = c[x]**2 * dt**2/ dx**2 *(-p[x-2]+16*p[x-1]-30*p[x]+16*p[x+1]-p[x+2])/12 +2*p[x] - p_old[x] 
    if boundary == 'zero':
        p_new[0:1] = 0
        p_new[-2:] = 0
    elif boundary == 'neumann':
        p_new[0:1] = p_new[2]
        p_new[-2:] = p_new[-3]
    # elif boundary == 'absorbing':
    #     p_new[0:1] = p_new[2:]
    #     p_new[-2:] = 0
    p_old = p.copy()
    p=p_new.copy()
    line1.set_ydata(p)
    ax.set_ylim(np.min(p), np.max(p))
    ax.set_title(f'CFL= {cfl:.2f}, Time = {t*dt:.3f} s')
    t+=1
    if t%5==0:
        return(line1,)
ani = animation.FuncAnimation(fig, func= five_pt_stencil, frames=np.arange(0,10000,5), interval=1, blit=False)
# Writer = writers['ffmpeg']
# writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)
# ani.save('1D_fd.mp4', writer=writer)
plt.show()