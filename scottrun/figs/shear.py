import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Affine2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
#sformatter=ScalarFormatter(useOffset=True,useMathText=False)
#sformatter.set_scientific(True)
#sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,7))
fig = plt.figure(1)
ax = fig.add_axes([0.16,0.12,0.8,0.43])
dPde=0.152181
Pressure=0.0404889
epsilon=0.24014

mydata = np.loadtxt('sheartest_rhoB0.dat',skiprows=1,unpack=True)
pizzPhydro=-mydata[0]
pioverPhydro=-mydata[1]
pizzcheck=mydata[2]
epsilon=mydata[3]

mydata2 = np.loadtxt('sheartest_rhoB0.1.dat',skiprows=1,unpack=True)
pizzPhydro2=-mydata2[0]
pioverPhydro2=-mydata2[1]
pizzcheck2=mydata2[2]
epsilon2=mydata2[3]

mydata3 = np.loadtxt('sheartest_deltaf_rhoB0.dat',skiprows=1,unpack=True)
pizzPhydro3=-mydata3[0]
pioverPhydro3=-mydata3[1]
pizzcheck3=mydata3[2]
epsilon3=mydata3[3]

mydata4 = np.loadtxt('sheartest_deltaf_rhoB0.1.dat',skiprows=1,unpack=True)
pizzPhydro4=-mydata4[0]
pioverPhydro4=-mydata4[1]
pizzcheck4=mydata4[2]
epsilon4=mydata4[3]

plt.plot(pioverPhydro,pizzcheck,marker='o',markersize='8',linestyle='--',linewidth=1,color='g')
plt.plot(pioverPhydro2,pizzcheck2,marker='s',markersize='8',linestyle='--',linewidth=1,color='r')
plt.plot(pioverPhydro3,pizzcheck3,marker='o',markersize='8',linestyle='--',linewidth=1,color='b')
plt.plot(pioverPhydro4,pizzcheck4,marker='s',markersize='8',linestyle='--',linewidth=1,color='k')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,1.5,0.2), minor=False)
ax.set_xticklabels(np.arange(0,1.5,0.2), minor=False)
ax.set_xticks(np.arange(0,1.5,0.05), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,1.0)

ax.tick_params(direction='in', which='both')
ax.yaxis.set_ticks_position('both')
ax.set_yticks(np.arange(0.5,1.5,0.1), minor=False)
ax.set_yticklabels(np.arange(0.5,1.5,0.1), minor=False)
ax.set_yticks(np.arange(0.5,1.5,0.05), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.ylim(0.58,1.05)

text(0.98,0.86,'$\\rho_B=0$',color='green',fontsize='22',ha='right')
text(0.98,0.76,'$\\rho_B=0.1$',color='red',fontsize='22',ha='right')
text(0.5,0.8,'$\\rho_B=0$',color='b',fontsize='22',ha='right')
text(0.78,0.66,'$\\rho_B=0.1$',color='k',fontsize='22',ha='right')

plt.xlabel('$-\pi^{\\rm(hydro)}_{zz}/P$', fontsize=18, weight='normal')
plt.ylabel('$\pi_{zz}/\pi^{\\rm(hydro)}_{zz}$',fontsize=18)


ax = fig.add_axes([0.16,0.55,0.8,0.43])
plt.plot(pioverPhydro,epsilon,marker='o',markersize='8',linestyle='--',linewidth=1,color='g')
plt.plot(pioverPhydro2,epsilon2,marker='s',markersize='8',linestyle='--',linewidth=1,color='r')
unity=epsilon2
unity=epsilon2+1.0-epsilon2
unityk=epsilon
unityk=epsilon+1.0-epsilon
plt.plot(pioverPhydro2,unity,marker='o',markersize='8',linestyle='--',linewidth=1,color='b')
plt.plot(pioverPhydro,unityk,marker='s',markersize='8',linestyle='--',linewidth=1,color='k')

ax.set_xticks(np.arange(0,1.5,0.2), minor=False)
ax.set_xticks(np.arange(0,1.5,0.05), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.set_xticklabels([], minor=False)
plt.xlim(0.0,1.25)

ax.tick_params(direction='in', which='both')
ax.yaxis.set_ticks_position('both')
ax.set_yticks(np.arange(0.4,1.2,0.02), minor=False)
ax.set_yticklabels(np.arange(0.4,1.2,0.02), minor=False)
ax.set_yticks(np.arange(0.4,1.2,0.01), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.ylim(0.99,1.07)
plt.ylabel('$\epsilon/\epsilon_{\\rm hydro}$',fontsize=18)

text(0.6,1.04,'$p_i=p^\prime_i-\\alpha \pi_{ij}p^\prime_j$',fontsize=20,ha='left')
text(0.74,1.005,'$\delta f/f\sim p_i\pi_{ij}p_j/E$',fontsize=20,ha='left')

plt.savefig('shear.pdf',format='pdf')
os.system('open -a Preview shear.pdf')
#plt.show()
quit()
