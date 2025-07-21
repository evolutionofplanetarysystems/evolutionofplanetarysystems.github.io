import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

import plotsetup

mass_min=.1
mass_max=4106
mass=np.logspace(np.log10(mass_min),np.log10(mass_max),1000)

# https://web.gps.caltech.edu/classes/ge131/notes2023/Chapter_9.pdf

m1=np.array([1,10,100,1000])
rock=np.array([1,1.5,2.2,2.8])
ice=np.array([1.2,2.5,4.2,4.8])
gas=np.array([2.6,6.1,9.5,11])

colors=plt.rcParams["axes.prop_cycle"].by_key()['color']


color_rock="brown"
color_ice="blue"
color_icegiant="green"
color_gas="black"

color_rock=color_ice=color_icegiant=color_gas="black"

fig=plt.figure(figsize=(8,5))
do_dots=False
if do_dots:
	plt.plot(m1,m1*5.97e27/(4*np.pi/3*(rock*6378e5)**3),'-o',color=color_rock)
	plt.plot(m1,m1*5.97e27//(4*np.pi/3*(ice*6378e5)**3),'-o',color=color_ice)
	plt.plot(m1,m1*5.97e27//(4*np.pi/3*(gas*6378e5)**3),'-o',color=color_gas)
plt.xscale("log");plt.yscale("log")

x=np.log10(mass)
yrock=.4*x +.2*(x/2)**2
rrock = 5.5*10**yrock
yice=.18*x**2-.18*x+.5
rice=10**yice
rice=np.clip(rrock*.2,a_min=1,a_max=np.inf)
ygas=.21*x**2-.24*x-.57
rgas=10**ygas
rgas[0:np.argmin(rgas)]=np.min(rgas)
rgas=np.clip(rgas,.2,np.inf)
# using exoplex:
a=.1 ; rho1=10 ; exponent=0.35
rmetal=rho1*((a/mass+1)/(1.+a)*mass)**exponent
a=.3 ; rho1=4.6 ; exponent=0.22
rrock=rho1*((a/mass+1)/(1.+a)*mass)**exponent
a=.1 ; rho1=2.1 ; exponent=0.2
rice=rho1*((a/mass+1)/(1.+a)*mass)**exponent
a=.4 ; rho1=1 ; exponent=0.25
rwater=rho1*((a/mass+1)/(1.+a)*mass)**exponent
radius_gg=(0.96 + 0.21*np.log10(mass*5.97/1900)-0.20*(np.log10(mass*5.97/1900))**2)*69911/6378
rgg=mass*5.97e27/(4*np.pi/3*(radius_gg*6378e5)**3)
l=np.where(mass<100)
plt.plot(mass[l],rmetal[l],linestyle="dotted",color="black")#,label="iron")
plt.plot(mass[l],rrock[l],color=color_rock)#,label="rock")
plt.plot(mass[l],rice[l],'--',color=color_ice)#,label="ice")
l=np.where(mass>8)
plt.plot(mass[l],rgas[l],linestyle="dashdot",color=color_gas)#,label="gas")

plt.xlabel("Planet Mass (Earth masses)")
plt.ylabel(r"Bulk Density (g/cm$^3$)")

namep=['M','V','E','J','S','U','N']
mp=np.array([0.64,4.87,5.97,1900,568,87,102])/5.97
rhop=np.array([3.93,5.2,5.5,1.33,0.69,1.32,1.64])
#plt.scatter(mp,rhop,color="black")
plt.scatter(mp,rhop,s=150,edgecolors="black",facecolors="white",zorder=5)
for n,m,rho in zip(namep,mp,rhop):
	if m>mass_min:
		plt.text(m,rho,n,va="center",ha="center",zorder=6,fontweight='demibold',fontsize=10)

#plt.plot(mass,4*mass**(1/2),color="grey",linestyle="dashed")
ricerock=(rice*rrock)**.5

rho_min=.1
alpha_fill=0 #.05
plt.fill(np.append(mass,mass[-1::-1]),np.append(1.5*rrock,rice[-1::-1]),color=color_rock,alpha=alpha_fill)
#plt.fill(np.append(mass,mass[-1::-1]),np.append(ricerock,rice[-1::-1]),color=color_ice,alpha=alpha_fill)
plt.fill(np.append(mass,mass[-1::-1]),np.append(rice,1.5*rgas[-1::-1]),color=color_icegiant,alpha=alpha_fill)
plt.fill(np.append(mass,[mass[-1],mass[0]]),np.append(1.5*rgas,[rho_min,rho_min]),color=color_gas,alpha=alpha_fill)

rot=12
plt.text(80,50,"METAL",color="black",rotation=rot,fontstyle='italic')
plt.text(80,18,"ROCK & METAL",color="black",rotation=rot,fontstyle='italic')
plt.text(80,7,"ROCK & ICE",color="black",rotation=rot,fontstyle='italic')
plt.text(80,1.6,"ICE & GAS",color="black",rotation=rot,fontstyle='italic')
plt.text(50,.2,"GAS ",color="black",rotation=rot,fontstyle='italic')
plt.xlim(mass.min(),mass.max())
plt.ylim(rho_min,100)

do_exoplan_nearby=True
if do_exoplan_nearby:
	d=ascii.read("exoplan.dat")
	d=d[d['Density']>0]
	# fix proxima cen b and 51 peg b
	dp=d[d['planet']=="prox_cen_b"]
	d=d[d['planet']!="prox_cen_b"]
	d=d[d['planet']!="51_peg_b"]
	plt.scatter(d['Mass'],d['Density'],20,label="Nearby exoplanets",color="black")
	plt.errorbar([1.1],[2],yerr=[1.5],lolims=True,fmt='o',capsize=2,color="black",markersize=5)
	plt.errorbar([0.5*1900/5.97],[0.2],yerr=[.1],lolims=True,fmt='o',capsize=2,color="black",markersize=5)
# EXOPLANETS
do_exoplan=True
if do_exoplan:
	de=ascii.read("exoplanet.eu_catalog_07-07-25_00_31_54.csv")
	de=de[de['orbital_period']>1]
	massp=de['mass']*1900e27
	mearth=massp/5.97e27
	radius=de['radius']*69911e5
	density = massp/(4*np.pi/3*radius**3)
	de['density']=density
	plt.scatter(mearth,density,1,label="Exoplanets",zorder=0,color="grey")
	dcrit=.1*(mearth/100)**(.95)*3

do_browndwarf=True
if do_browndwarf:
	d=ascii.read("nc+0.0_co1.0_mass")
	m=d['M/Msun']
	r=d['R/Rsun']
	density=m*1.99e33/(4*np.pi/3*(r*6.96e10)**3)
	t=d['age(Gyr)']
	mme=m*1.99e33/5.97e27
	l=np.where(t==4.0)
	#plt.plot(mme[l],density[l]*1.1,color="grey")
	l=np.where(t==0.01)
	plt.plot(mme[l],density[l],linestyle="dashdot",color="grey")


#plt.legend()
plt.tight_layout()
plt.savefig("planmod.pdf")

print("WHAT TO DO WITH SIZE...DONT CALL SMALL PLANETS GIANTS")
