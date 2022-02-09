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
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(5,10))
fig = plt.figure(1)

mydata = np.loadtxt('rescount_polemass.dat',skiprows=0,unpack=True)
rescount_pm=mydata[4]
restarget_pm=mydata[3]
resratio_pm=mydata[2]
resmass_pm=mydata[1]
error_pm=1.0/sqrt(rescount_pm)

mydata = np.loadtxt('rescount_SF.dat',skiprows=0,unpack=True)
rescount_sf=mydata[4]
restarget_sf=mydata[3]
resratio_sf=mydata[2]
resmass_sf=mydata[1]
error_sf=1.0/sqrt(rescount_sf)

ax = fig.add_axes([0.19,0.15,0.77,0.27])

plt.plot(resmass_pm,resratio_sf,linestyle='None',color='r',marker='o',markersize=2)
plt.errorbar(resmass_pm,resratio_sf,yerr=error_pm,linestyle='None',color='r',linewidth=1)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,4.0,0.5), minor=False)
ax.set_xticklabels(np.arange(0,4.0,0.5), minor=False, family='serif')
ax.set_xticks(np.arange(0,4.0,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0,2.5)

ax.set_yticks(np.arange(0.8,1.2,0.05), minor=False)
ax.set_yticklabels(np.arange(0.8,1.2,0.05), minor=False, family='serif')
ax.set_yticks(np.arange(0.8,1.2,0.01), minor=True)
plt.ylim(0.9,1.125)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
#ax.yaxis.set_major_formatter(sformatter)
plt.xlabel('$M$ (GeV/$c^2$)', fontsize=18, weight='normal')
plt.ylabel('$N_{\\rm sampler}/N_{\\rm therm}$',fontsize=18)

text(0.1,0.95,'pole mass',color='r',size=20)

#############################################

ax = fig.add_axes([0.19,0.42,0.77,0.27])

plt.plot(resmass_sf,resratio_pm,linestyle='None',color='b',marker='o',markersize=2)
plt.errorbar(resmass_sf,resratio_pm,yerr=error_sf,linestyle='None',color='b',linewidth=1)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,4.0,0.5), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,4.0,0.1), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0,2.5)

ax.set_yticks(np.arange(0.8,1.2,0.05), minor=False)
ax.set_yticklabels(np.arange(0.8,1.2,0.05), minor=False, family='serif')
ax.set_yticks(np.arange(0.8,1.2,0.01), minor=True)
plt.ylim(0.9,1.125)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
#ax.yaxis.set_major_formatter(sformatter)
plt.xlabel(None)
plt.ylabel('$N_{\\rm sampler}/N_{\\rm therm}$',fontsize=18)

text(0.1,0.95,'spectral function',color='b',size=20)

#############################################

ax = fig.add_axes([0.19,0.69,0.77,0.27])

#error_ratio=sqrt(error_sf*error_sf+error_pm*error_pm)
resratio=restarget_sf/restarget_pm

plt.plot(resmass_sf,resratio,linestyle='None',color='g',marker='o',markersize=2)
#plt.errorbar(resmass_sf,resratio,yerr=error_ratio,linestyle='None',color='g',linewidth=1)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,4.0,0.5), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,4.0,0.1), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0,2.5)

ax.set_yticks(np.arange(0.0,3,0.5), minor=False)
ax.set_yticklabels(np.arange(0.0,3,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(0.0,3,0.1), minor=True)
plt.ylim(0.25,2.5)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
#ax.yaxis.set_major_formatter(sformatter)
plt.xlabel(None)
plt.ylabel('$N_{\\rm SF}/N_{\\rm PM}$',fontsize=18)

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('rescount.pdf',format='pdf')
os.system('open -a Preview rescount.pdf')
#plt.show()
quit()
