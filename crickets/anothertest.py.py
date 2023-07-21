# import numpy as np
# import matplotlib.pyplot as plt

# x = np.arange(0,25)
# y = np.sin(x)

# plt.plot(x,y)

# plt.savefig(fname='/home/sofairj/test_plot.pdf', dpi=300, format='pdf')
# plt.savefig(fname='/home/sofairj/test_plot.png', dpi=300, format='png')
# plt.savefig(fname='/home/sofairj/test_plot.jpg', dpi=300, format='jpg')


# a = 'This is a string with/ repetitive characters,/ especially forward/ slashes/'
# b='////'

# print(b.rfind('a'))

# import os

# cfrk = os.getenv("cfrk")
# # print(f'\n\n{cfrk}', type(cfrk))
# print(f"\n\n{os.path.realpath(os.path.expanduser(f'~/{cfrk}'))}")
# # print(f"\n\n{os.path.expanduser('~/nrao_summer_2023')}")

from nrao_summer_2023.cricket.cricket.plotting import save_fig
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0,10)
y = np.sin(x)

fig, ax = plt.subplots()
ax.plot(x,y)

save_fig('this is a file', types=['png', 'pdf'])