import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import subprocess
# import mayavi.mlab as mlab
# import time as timel
### Initial conditions
nx= int(250) ## Number of x grid points
ny= int(250) ## Number of y grid points
nt= 1000 ## Number of time steps
dx=1
dy=1
dt=0.001
time2= np.linspace(0,nt*dt,nt+1)
p_old,p,p_new = (np.zeros((nx,ny)) for _ in range(3))
c0=700
c = np.zeros((nx,ny)) + c0
domain_x,domain_y =np.meshgrid(np.arange(0, nx*dx, dx), np.arange(0, ny*dy, dy))
### Source time function
f0=100.
t0=4./f0
time = np.linspace(0*dt, nt*dt, nt)

def src_func(time):
    return -8. * (time - t0) * f0 * (np.exp(-1.0 * (4*f0) ** 2 * (time - t0) ** 2))
isrc = [round(nx/2), round(ny/2)]  # Source location as a list for 2D
### Boundary conditions
boundary = 'zero'  

### Plot setup
cfl = c0 * dt * np.sqrt(1/dx**2 + 1/dy**2)
if cfl > 1:
    print(f'Warning: CFL condition violated! CFL = {cfl:.2f} > 1.0')
else:
    print(f'CFL condition satisfied. CFL = {cfl:.2f} <= 1.0')
# fig,ax = plt.subplots(figsize=(10, 5))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')# line1= ax.imshow(p)
plot = [ax.plot_surface(domain_x, domain_y, p, cmap='magma', edgecolor='none')]
# ### Solver
t=0
z =np.zeros((nx, ny,nt))
# def five_pt_stencil(plot):
#     global p,p_old,p_new, t
for frame in np.arange(0,nt): 
    for x in np.arange(2, nx-2, 1):
        for y in np.arange(2, ny-2, 1):
            if x == isrc[0] and y == isrc[1]:
                p_new[x,y] = c[x,y]**2 * (dt**2/ dx**2 *(-p[x-2,y]+16*p[x-1,y]-30*p[x,y]+16*p[x+1,y]-p[x+2,y])/12 +dt**2/ dy**2 *(-p[x,y-2]+16*p[x,y-1]-30*p[x,y]+16*p[x,y+1]-p[x,y+2])/12) + (1e-1 *src_func(t*dt)/dx * dt**2) +2*p[x,y] - p_old[x,y] 
            else:
                p_new[x,y] = c[x,y]**2 * (dt**2/ dx**2 *(-p[x-2,y]+16*p[x-1,y]-30*p[x,y]+16*p[x+1,y]-p[x+2,y])/12 +dt**2/ dy**2 *(-p[x,y-2]+16*p[x,y-1]-30*p[x,y]+16*p[x,y+1]-p[x,y+2])/12) +2*p[x,y] - p_old[x,y]
        if boundary == 'zero':
            p_new[0:1,:] = 0
            p_new[-2:,] = 0
            p_new[:,0:1] = 0
            p_new[:,-2:] = 0
        elif boundary == 'neumann':
            p_new[0,:] = p_new[1,:]
            p_new[-1,:] = p_new[-2,:]
            p_new[:,0] = p_new[:,1]
            p_new[:,-1] = p_new[:,-2]
    p_old = p.copy()
    p=p_new.copy()
    t += 1
    # plot[0].remove()
    # plot[0] = ax.plot_surface(domain_x, domain_y, p, cmap='magma', edgecolor='none')
    # ax.set_title(f'Time = {t*dt:.3f} s')
    z[:,:,frame] = np.copy(p)

def update_plot(frame, z):
    ax.clear()
    ax.set_zlim(-0.05, 0.05)  # Optional: set consistent z-axis limits
    surf = ax.plot_surface(domain_x, domain_y, z[:, :, frame], cmap='magma', edgecolor='none')
    ax.set_title(f'Time = {frame*dt:.3f} s')
    return [surf]
# ani = animation.FuncAnimation(fig, func= five_pt_stencil, frames=np.arange(0,1000,5), interval=10, blit=False)
# fn = 'surface_2D_fd'
# ani.save(fn+'.mp4',writer='ffmpeg',fps=20)
# ani.save(fn+'.gif',writer='imagemagick',fps=20)
# cmd = 'magick convert %s.gif -fuzz 5%% -layers Optimize %s_r.gif'%(fn,fn)
# subprocess.check_output(cmd)

ani = animation.FuncAnimation(
    fig, update_plot, frames=np.arange(0, 1000, 5),
    fargs=(z,), interval=1, blit=False
)
plt.show()
       