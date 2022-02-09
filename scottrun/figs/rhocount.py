import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.17,0.14,0.80,0.85])

mydata = np.loadtxt('rhocount_polemass.dat',skiprows=0,unpack=True)
T=mydata[0]
rhocount=mydata[1]
deltacount=mydata[2]
mydataSF = np.loadtxt('rhocount_SF.dat',skiprows=0,unpack=True)
TSF=mydataSF[0]
rhocountSF=mydataSF[1]
deltacountSF=mydataSF[2]
#yscale('log')
plt.plot(T,rhocount,linestyle='-',linewidth=2,color='r')
plt.plot(T,deltacount,linestyle='-',linewidth=2,color='r')
plt.plot(TSF,rhocountSF,linestyle='-',linewidth=2,color='g')
plt.plot(TSF,deltacountSF,linestyle='-',linewidth=2,color='g')

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(80,200,20), minor=False)
ax.set_xticklabels(np.arange(80,200,20), minor=False, family='serif',fontsize=18)
ax.set_xticks(np.arange(80,200,10), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.xlim(99,161)

ax.set_yticks(np.arange(0.0,0.05,0.01), minor=False)
ax.set_yticklabels(np.arange(0,0.05,0.01), minor=False, family='serif',fontsize=18)
ax.set_yticks(np.arange(0,0.05,0.005), minor=True)
plt.ylim(0.0,0.0375)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
#ax.yaxis.set_major_formatter(sformatter)

text(141,0.0005,'$\Delta$ baryons',fontsize='22',ha='left')
text(142,0.0135,'$\\rho$ mesons',fontsize='22',ha='left')

text(100,0.027,'Using pole masses',fontsize='22',color='red',ha='left')
text(100,0.024,'SMASH Spectral Functions',fontsize='20',color='green',ha='left')

plt.xlabel('$T$ (MeV)', fontsize=22, weight='normal')
plt.ylabel('$\langle n\\rangle$ (fm$^{-3}$)',fontsize=24,labelpad=-2)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('rhocount.pdf',format='pdf')
os.system('open -a Preview rhocount.pdf')
#plt.show()
quit()
