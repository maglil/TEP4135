import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import pi

fig =  plt.figure()
ax = plt.axes(xlim = [0,2*pi],ylim = [-1,1]) #axes is the key thing were things go. has helper methods

#Defin

ln, = ax.plot([],[],'ko')

def init():
    x = 0
    y = 0
    ln.set_data(x,y)
    #ln, = ax.plot(x,y,'ko')
    return ln,

def update(i):
    x = 2*pi*i/50
    y = np.sin(x)
    ln.set_data(x,y)
    return ln,

ani = FuncAnimation(fig, update, frames=50, interval=20, init_func=init, blit=True)

plt.show()
