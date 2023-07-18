import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0,25)
y = np.sin(x)

plt.plot(x,y)

plt.savefig(fname='/home/sofairj/test_plot.pdf', dpi=300, format='pdf')
plt.savefig(fname='/home/sofairj/test_plot.png', dpi=300, format='png')
plt.savefig(fname='/home/sofairj/test_plot.jpg', dpi=300, format='jpg')