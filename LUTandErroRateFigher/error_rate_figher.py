import math
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)

np.set_printoptions(precision=15)
q = 2 ** 32
q_2 = q*q
B = 4
B8 = 8
wBML1 = 2 * 1024/(2**1)
wBML2 = 2 * 1024/(2**2)
wBML3 = 2 * 1024/(2**3)

wFB = 2 * 1024
n = 630
a4 = n/48
y1 = np.empty(0)
y2 = np.empty(0)
y3 = np.empty(0)
y4 = np.empty(0)

y1_8 = np.empty(0)
y2_8 = np.empty(0)
y3_8 = np.empty(0)
y4_8 = np.empty(0)

xBML = np.empty(0)
xBML2 = np.empty(0)
xBML3 = np.empty(0)

xBML_8 = np.empty(0)
xBML2_8 = np.empty(0)
xBML3_8 = np.empty(0)

ylevel = np.empty(0)
xlevel = np.empty(0)

xFB = np.empty(0)
xFB_8 = np.empty(0)

c2 = q_2/(12*wFB*wFB)
c3 = n*q_2/(24*wFB*wFB)
for i in range( 3000,8000):
    taro = i/1000
    Proby = math.erf(taro/(2**(1/2)))
    a1 = q_2/(16*B*B*taro*taro)
    a1_8 = q_2/(16*B8*B8*taro*taro)

    bBML_x = a1 - q_2/(12*wBML1*wBML1) + 1/12 - n*q_2/(24*wBML1*wBML1) -a4
    bBML_x2 = a1 - q_2/(12*wBML2*wBML2) + 1/12 - n*q_2/(24*wBML2*wBML2) -a4
    bBML_x3 = a1 - q_2/(12*wBML3*wBML3) + 1/12 - n*q_2/(24*wBML3*wBML3) -a4

    bBML_x_8 = a1_8 - q_2/(12*wBML1*wBML1) + 1/12 - n*q_2/(24*wBML1*wBML1) -a4
    bBML_x2_8 = a1_8 - q_2/(12*wBML2*wBML2) + 1/12 - n*q_2/(24*wBML2*wBML2) -a4
    bBML_x3_8 = a1_8 - q_2/(12*wBML3*wBML3) + 1/12 - n*q_2/(24*wBML3*wBML3) -a4

    bFB_x = a1 - c2 + 1/12 - c3 -a4
    bFB_x_8 = a1_8 - c2 + 1/12 - c3 -a4


    if (bBML_x > 0):
        xBML = np.append(xBML, bBML_x)
        y1 = np.append(y1, Proby)
    if (bBML_x2 > 0):
        xBML2 = np.append(xBML2, bBML_x2)
        y2 = np.append(y2, Proby)
    if (bBML_x3 > 0):
        xBML3 = np.append(xBML3, bBML_x3)
        y3 = np.append(y3, Proby)
    if (bFB_x > 0):
        xFB = np.append(xFB, bFB_x)
        y4 = np.append(y4, Proby)

    if (bBML_x_8 > 0):
        xBML_8 = np.append(xBML_8, bBML_x_8)
        y1_8 = np.append(y1_8, Proby)
    if (bBML_x2_8 > 0):
        xBML2_8 = np.append(xBML2_8, bBML_x2_8)
        y2_8 = np.append(y2_8, Proby)
    if (bBML_x3_8 > 0):
        xBML3_8 = np.append(xBML3_8, bBML_x3_8)
        y3_8 = np.append(y3_8, Proby)
    if (bFB_x_8 > 0):
        xFB_8 = np.append(xFB_8, bFB_x_8)
        y4_8 = np.append(y4_8, Proby)

    xlevel = np.append(xlevel, 2**17)
    ylevel = np.append(ylevel, Proby)


plt.figure(dpi=100)
plt.plot(xFB, y4,  c='k',linestyle='-', label='TreeFB')
plt.plot(  xBML,y1, c='orange',linestyle='-', label='$\\vartheta$ = 1')
plt.plot(  xBML2,y2, c='lawngreen',linestyle='-', label='$ \\vartheta$ = 2')
plt.plot(  xBML3,y3, c='cyan',linestyle='-', label='$ \\vartheta$ = 3')
plt.plot(  xlevel,ylevel, c='red',linestyle='--', label='127-bit lwe Security')

plt.plot(xFB_8, y4_8, c='k',linestyle='-.')
plt.plot(  xBML_8,y1_8, c='orange',linestyle='-.')
plt.plot(  xBML2_8,y2_8, c='lawngreen',linestyle='-.')
plt.plot(  xBML3_8,y3_8, c='cyan',linestyle='-.')
plt.plot(  [],[], c='k',linestyle='-.', label='denote B = 8')


new_ticks = [7*pow(10,15), 6*pow(10,15), 5*pow(10,15), 4*pow(10,15), 3*pow(10,15), 2*pow(10,15), 1*pow(10,15), 0]
plt.xticks(new_ticks, ['7', '6', '5', '4', '3', '2', '1', '0'])

#plt.xlim(7*pow(10,15),-1000)



plt.text(xlevel[0]-0.2*10**15, ylevel[0],'$2^{-15}\cdot q$', va = 'top')


#ax=plt.axes()
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)


#print(xFB)
#print(xBML2)
#print(xBML3)

plt.xlabel('error varience $ (10^{15})$', fontsize=14)
plt.ylabel('correctness of the ModSwitch', fontsize=14)
plt.legend()

plt.show()
