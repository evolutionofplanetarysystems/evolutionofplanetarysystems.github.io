from astropy.io import ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import plotsetup

mods=Table()
mods['file']=['0005000M.track.eep','0010000M.track.eep','0020000M.track.eep']
mods['label']=["0.5","1","2"]
mods['age_max']=[13e9,11.45e9,1.35e9]
mods['Tms']=[3900,5900,9200]
mods['Lms']=[.02,.7,12]
mods['tms1']=[0.23,0.1,0.1]
mods['tms2']=[13,10.,1.06]
mods['trg1']=[13,11.34,1.15]
mods['trg2']=[0,11.44,1.32] # includes AGB

# solsys now
rhab1=0.95 ;  rhab2=1.67

plt.clf()
plt.xlim(0.23,12.8)
plt.ylim(0,5)

for imod,mod in enumerate(mods):
	d=ascii.read(mod['file'])
	age=d['star_age']/1e9
	Teff=10**d['log_Teff']
	L=10**d['log_L']
	if imod==1: L=L/1.102
	R=10**d['log_R']
	l=np.where((age>mods[imod]['tms1'])&(age<mods[imod]['tms2']))
	age=age[l] ; L=L[l]
	rhab11=rhab1*L**.25
	rhab22=rhab2*L**.25
	plt.plot(age,rhab11,color="black")
	plt.plot(age,rhab22,color="black")
	plt.fill(np.append(age,age[-1::-1]),np.append(rhab11,rhab22[-1::-1]),color="grey",alpha=.2)
plt.xlabel("age Gyr")
plt.ylabel("Distance from star (au)")
plt.scatter([4.5],[1],color="black")
plt.text(11,.58,r"$0.5 M_\odot$",ha="center",va="center")
plt.text( 8,1.4,r"$1 M_\odot$",ha="center",va="center")
plt.text(.7,2.7,r"$2 M_\odot$",ha="center",va="center")
plt.tight_layout()
plt.savefig("msvary.pdf")