from astropy.io import ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


mods=Table()
mods['file']=['track0p5.csv','track1p0.csv','track2p0.csv']
mods['label']=["0.5","1","2"]
mods['age_max']=[13e9,11.45e9,1.35e9]
mods['Tms']=[3900,5900,9200]
mods['Lms']=[.02,.7,12]
mods['tms1']=[0.1,0.1,0.1]
mods['tms2']=[13,10.,1.06]
mods['trg1']=[13,11.34,1.15]
mods['trg2']=[0,11.44,1.32] # includes AGB

# diagnostic plots to check if the boundaries between long-lived phases are correct
for mod in mods:
	d=ascii.read(mod['file'])
	age=d['star_age']/1e9
	Teff=10**d['log_Teff']
	L=10**d['log_L']
	R=10**d['log_R']
	fig,ax=plt.subplots(3)
	ax[0].plot(age,Teff) ; ax[0].set_ylabel('Teff')
	ax[1].plot(age,L)    ; ax[1].set_ylabel('L')
	ax[2].plot(age,R)    ; ax[2].set_ylabel('R')
	ax[2].set_xlabel("age Gyr")
	ax[0].set_ylim(3000,16000)
	ax[2].set_ylim(0,30)
	for aax in ax:
		aax.axvline(mod['tms1'],linestyle='dotted',color='black')
		aax.axvline(mod['tms2'],linestyle='dotted',color='black')
		aax.axvline(mod['trg1'],linestyle='dotted',color='black')
		aax.axvline(mod['trg2'],linestyle='dotted',color='black')
		aax.set_xlim(.1,age.max())
		#aax.set_xlim(9,11.5)
	plt.tight_layout()

fig = plt.figure(figsize=(6, 6)) 
plt.rcParams['axes.labelsize'] = 12
ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', labelsize=10)
plt.rc('font',  family='serif', size=12)
plt.rc('text',  usetex=False)
for spine in ax.spines.values():
        spine.set_linewidth(2)  # Set the axis thickness
colors=plt.rcParams["axes.prop_cycle"].by_key()['color']
import matplotlib.ticker as ticker
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}')) # with commas

age_min=0.01

for i,mod in enumerate(mods):
	d=ascii.read(mod['file'])
	d=d[(d['star_age']<mod['age_max'])&(d['star_age']>age_min)]
	age=d['star_age']/1e9
	Teff=10**d['log_Teff']
	L=10**d['log_L']
	R=10**d['log_R']
	print(mod['label'],age.min(),age.max())
	plt.plot(Teff,L,label=mod['label'],color=colors[i],linewidth=1)
# phases
	for t1,t2,phase in zip([mod['tms1'],mod['trg1']],[mod['tms2'],mod['trg2']],['main_seq','red_giant']):
		l=np.where((age>t1)&(age<t2))
		if len(l[0])>0:
			plt.plot(Teff[l],L[l],color=colors[i],linewidth=4)
			print(f"{phase:12}: {t1:5.2f}-{t2:5.2f} Gyr, {int(Teff[l].min()):6}-{int(Teff[l].max()):6} K, "+
				f"{R[l].min():6.2f}-{R[l].max():6.2f} Rs  {L[l].min():6.1f}-{L[l].max():6.1f} Ls")
		else:
			print(f"{phase} does not occur, age range {t1:5.2f}-{t2:5.2f} Gyr")	
# start and end
	plt.text(Teff[0],L[0],r'$\star$',
		color=colors[i],ha="center",va="center",fontsize=30)
	plt.text(Teff[-1],L[-1],r'$\circ$',
		color=colors[i],ha="center",va="center",fontsize=30)
	do_mslab=True
	if do_mslab:
		plt.text(mod['Tms'],mod['Lms'],mod['label']+r'$M_\odot$',
			color=colors[i],ha='right',va='center')
	do_agelabs=False
	if do_agelabs:
		agelabs=[1,2,4,8]
		agelabs=[mod['age_ms']]
		for agelab in agelabs:
			if agelab<age.max():
				l=np.where(age<agelab)
				ilab=l[0][-1]
				plt.text(10**d[ilab]['log_Teff'],10**d[ilab]['log_L'],str(agelab),
					color=colors[i],ha="center",va="center",fontsize=20)
# ZAMS
d=ascii.read("isochrone_8.csv")
d=d[d['initial_mass']<4]
Tzams=10**d['log_Teff'] ; Lzams=10**d['log_L']
plt.plot(Tzams,Lzams,color="black")
plt.text(7600,4.7,'ZAMS',color="black",rotation=-28,ha='center',va='center')
plt.text(9200,80,"Main\nLifetime",color="black")
plt.fill(np.append(Tzams,Tzams[-1::-1]),np.append(Lzams,Lzams[-1::-1]*(1+Lzams[-1::-1]**1.2)),alpha=.2,color="grey")
# RGB
plt.show();ax.add_patch(patches.Ellipse((4600,450), 3200, 880, angle=0, facecolor='red', alpha=0.2))
plt.text(5600,90,"Red\nGiant",color="red")

plt.xlim(10000,3000)
plt.ylim(0.01,1000)
plt.yscale("log")
if not(do_mslab): plt.legend(title='Mass (solar)')
plt.xlabel(r"Temperature (K)")
plt.ylabel(r"Luminosity (solar)")

plt.tight_layout()
plt.savefig("tracks.pdf")