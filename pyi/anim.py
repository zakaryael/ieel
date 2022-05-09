from pylab import *
from matplotlib import animation
from utils import *

def init():
    ax.set_xlim(-3, 2)
    ax.set_ylim(-2.5, 2.5)
    del xdata[:]
    del ydata[:]
    line.set_data(xdata, ydata)
    head.set_data(xdata, ydata)
    return line,

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2, color='grey', alpha=0.7)
head, = ax.plot([], [], '.', lw=0.2, color='grey', alpha=0.9)
ax.grid()
xdata, ydata = [], []


def animate(i): 
    state = get_state(i, '../output_data/pg3')
    y = state[0][1]
    x = state[0][0]
    ax.set_xlim(x.mean()-3, x.mean() + 2)
    ymin, ymax = ax.get_ylim()
    if y.max() > ymax or y.min() < ymin:
        ax.set_ylim(y.mean() - 2.5, y.mean() + 2.5)
        ax.figure.canvas.draw()
    line.set_data(x, y)
    head.set_data(x[0], y[0])
    return line, head
 
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=range(0,398000,10), interval=10)

show()