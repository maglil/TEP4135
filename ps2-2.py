import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

x0 = -1
y0 = -1
p = np.array([x0,y0]) # position of particle

fig = plt.figure()
ax = plt.axes()
line, = ax.plot([], [], lw=3)

def init():
    line.set_data([], [])
    x0 = -1
    y0 = -1
    p = np.array([x0,y0]) # position of particle
    return line,

def animate(i):
    global p
    vx = np.cos(p[0])**2+5
    vy = np.cos(p[1])
    v = np.array([vy,vx])
    p = p + v
    
    line, = ax.plot(p[0],p[1],'o')
    return line,

anim = FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

#plt.show()

