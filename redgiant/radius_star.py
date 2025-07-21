from astropy.io import ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import scipy.integrate as integrate

mods=ascii.read("mods.dat")
mods=mods[1:3]
'''
mods=Table()
mods['file']=['track0p5','track1p0','track2p0']
mods['label']=["0.5","1","2"]
mods['age_max']=[13e9,11.45e9,1.353e9]
mods['Tms']=[3900,5900,9200]
mods['Lms']=[.02,.7,12]
mods['tms1']=[0.1,0.1,0.1]
mods['tms2']=[13,10.,1.06]
mods['trg1']=[13,11.34,1.15]
mods['trg2']=[0,11.44,1.32] # includes AGB
mods=mods[1:3]
'''

# diagnostic plots to check if the boundaries between long-lived phases are correct
fig,ax=plt.subplots(len(mods),figsize=(6,6))
letters='abcdef'

plt.rcParams['axes.labelsize'] = 12
for aax in ax: aax.tick_params(axis='both', labelsize=10)
plt.rc('font',  family='serif', size=12)
plt.rc('text',  usetex=False)
#for spine in ax.spines.values():
#        spine.set_linewidth(2)  # Set the axis thickness
colors=plt.rcParams["axes.prop_cycle"].by_key()['color']

color_R='black'
color_mdot=colors[3]

for i in range(len(mods)):
	mod=mods[i]
	d=ascii.read(mod['file']+'.csv')
	age=d['star_age']/1e9
	#agemin=mod['tms2']-0.1*(mod['tagb']-mod['tms2'])
	agemin=mod['tms2']
	d=d[(age<mod['age_max'])&(age>agemin)]
	age=d['star_age']/1e9
	Teff=10**d['log_Teff']
	L=10**d['log_L']
	R=10**d['log_R']
	mdot=-d['star_mdot']*1e6
	mass=d['star_mass']
	ax[i].plot(age,R,color=color_R)
	#ax[i].set_ylim(0.5,300)
	ax[i].set_yscale("log")
	ax[i].set_xlim(age.min(),age.max()+0.02*(mod['tagb']-mod['tms2']))
	lmax=np.argmax(mdot)
	p=interpolate.interp1d(age,mdot)
	mlost=integrate.quad(p,age.min(),age.max())[0]
	print(mod['file'], R.min(),R.max(),age[lmax],mdot[lmax],mlost)
	subcap='' #  f"({letter[i]})"
	ax[i].text(0.05,0.85,subcap+f"{mod['label']} "+r"$M_\odot$",transform=ax[i].transAxes)
	ax2 = ax[i].twinx()
	ax2.tick_params(axis='both', labelsize=10)
	#ax2.plot(age,mdot,color=color_mdot,alpha=.1)
	ax2.set_yscale("log")
	ax2.tick_params(axis='y', labelcolor=color_mdot)
	#ax2.set_ylim(mdot.max()*0.0007,mdot.max())
	mdotmin=1e-4
	ax2.set_ylim(mdotmin,5)
	ax2.fill(np.append(age,age[-1::-1]),np.append(mdot,mdot*0+mdotmin),color=color_mdot,alpha=.75)
	ax[i].set_ylim(R.min()-.1,R.max()+1)
	yl=ax[i].get_ylim() ; ycen=(yl[0]*yl[1])**.25
	xl=ax[i].get_xlim()
	for col,lab in zip(['tms2','trg','trc','tagb'],['Main','Red Giant ','','AGB ']):
		if mod[col]>xl[0]:
			ax[i].axvline(mod[col],color="grey",linestyle='dotted')
			ax[i].text(mod[col],ycen,lab,rotation=50,fontsize=8,va='bottom',ha='right')
ax[-1].set_xlabel("age (Gyr)")
#ax[-1].set_ylim(0,30)

plt.tight_layout()
aax=plt.gca()
fig.supylabel("       Radius (solar)")
#aax.text(-.1,.7,"Radius (solar radii)",rotation=90,transform=aax.transAxes)
aax.text(1.09,.4,"Mass loss rate (solar mass/Myr)",rotation=270,transform=aax.transAxes,color=color_mdot)

plt.savefig("radius_star.pdf")

yn=input("enter 'y' to get the table of maximum radii")
if yn=='y':
	plt.figure()
	import glob
	files=glob.glob('0*csv')
	files.sort()
	print(" M       dM         Rms      Tmax      Rmax     RmaxAU")
	Ms=[] ; Rmaxs=[] ; Tmaxs=[]
	for file in files:
		d=ascii.read(file)
		age=d['star_age']/1e9
		d=d[age>0.1]
		age=d['star_age']/1e9
		m=d['star_mass'] ; R=10**d['log_R']
		p=interpolate.interp1d(age,R)
		t=np.linspace(age.min(),age.max(),1000)
		RR=p(t)
		Rms=np.median(RR)
		ims=np.where(R>(1.5*R[0]))[0]
		lmax=np.argmax(R)

		l=np.where(R>215)
		if len(l[0])>0:
			age_engulf=age[l]
			print('1 au engulfment at ',age_engulf.min(),' lasting ',age_engulf.max()-age_engulf.min())

		print(round(m[0],2),' & ',round(m[0]-m[-1],4),' & ',round(Rms,3),' & ',round(age[lmax],2),' & ',round(R.max()),' & ',round(R.max()*6.96e10/1.496e13,2),' \\\\')
		#plt.plot(age,R)
		Ms=np.append(Ms,m[0])
		Rmaxs=np.append(Rmaxs,R.max())
		Tmaxs=np.append(Tmaxs,age[lmax])
		
	
	plt.scatter(Ms,Rmaxs)
	xx=np.linspace(.4,3,1000)
	plt.plot(xx,450/(1+28*np.exp(-4.5*xx)))
	plt.xlabel("Stellar Mass")
	plt.ylabel("Maximum radius (solar)")