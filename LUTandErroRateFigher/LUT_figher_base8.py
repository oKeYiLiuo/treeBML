import matplotlib.pyplot as plt
import numpy as np


x = [6,8,10,12]

y0 = [0.07132, 0.41267, 2.01964, 8.06182]
y1 = [0.04506, 0.20004, 0.89564, 3.46230]
y2 = [0.04898, 0.13760, 0.56824, 2.16922]
y3 = [0.05773, 0.13783, 0.41830, 1.51332]

plt.figure(dpi=100)
plt.plot(x, y0, marker='o', c='k',linestyle='-', label='TreeFB')
plt.plot(x, y1, marker='*' ,c='orange',linestyle='-.', label='Optimization with $\\vartheta$ = 1')
plt.plot(x, y2, marker='^' ,c='lawngreen',linestyle='--', label='Optimization with $\\vartheta$ = 2')
plt.plot(x, y3, marker='s' ,c='cyan',linestyle='dotted', label='Optimization with $\\vartheta$ = 3')

plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
new_ticks = [6,8,10,12]

plt.xticks(new_ticks, ['6-bit', '8-bit', '10-bit', '12-bit'])
plt.xlabel('LUT size')

plt.ylabel('running times/s ')

plt.legend()


plt.show()