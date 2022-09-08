import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig =  plt.figure()
ax = plt.axes() #axes is the key thing were things go. has helper methods
ln, = ax.plot([1,2],[2,1])


def init():    
    return ln,

def update(frame):    
    return ln,

ani = FuncAnimation(fig, update, frames=10, init_func=init, blit=True)

plt.show()
