from astropy.io import ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

d=ascii.read("0010000M.track.eep")

t=d['star_age']
mdot=-d['star_mdot']

plt.clf()
plt.plot(t,mdot,color="black")
plt.xlim(1.1462e10,1.1463e10)

g=Table()
g['amp']=np.array([9.3e-6,8.9e-6,1.5e-6,6e-7,5.2e-7,1.25e-5,1.7e-6,
	2.8e-5,1.3e-6,2.8e-6,
	1.9e-6
	])
g['cen']=np.array([1.14624100,1.14624115,1.1462413314,1.14623895,1.14626925,1.1462706435,1.14627083,
	1.14628332,1.146281,1.1462834815,
	1.146260018
	])*1e10
g['fwhm']=np.array([716,862,2000,65000,49000,715,1200,
	900,35000,1760,
	63000
	])
dmg=1.08*g['fwhm'].value*g['amp'].value

f=(8.*np.log(2))**.5
def mdeval(t,g,t0=1.146285e10,dt=89000,doplot=False):
	gsum=t*0
	for gg in g: 
		gi=gg['amp']*np.exp(-(t-gg['cen'])**2/2/(gg['fwhm']/f)**2)
		gsum=gsum+gi
		if doplot: plt.plot(t,gi)
	mdlast=.73e-15*np.clip(t-t0,a_min=0,a_max=dt)**2
	mdlast[(t-t0)>dt]=0
	troll=6000
	delt=(t-t0) - (dt-troll)
	l=np.where(delt>0)
	rf=t*0+1
	rf[l]=np.exp(-delt[l]**2/.2/troll**2)
	mdlast=mdlast*rf
	mdtot=gsum+mdlast
	if doplot: plt.plot(t,mdlast)
	return mdtot
mdtot=mdeval(t,g,doplot=True)
plt.plot(t,mdtot,color="gray")

print("Gaussian table")
ff=(8.*np.log(2))**.5
print("Amp(1e-6 Msolar/yr)  Cen(Myr)   Sigma (Myr)")
for i,gg in enumerate(g):
	print(f"{i:<2} & {gg['amp']*1e6:5.2f} & {(gg['cen']-11.462e9)/1e6:9.6f} & {gg['fwhm']/ff/1e6:9.6f} \\\\")
print("where cen is relative to time {tstart} yr")

print("arrays for fortran:")
print(g['amp'].value*1e6)
print((g['cen'].value-11.462e9)/1e6)
print(g['fwhm'].value/ff/1e6)

print(f"integrated lass in track = {np.trapz(mdot,t/1e10)*1e10}")
print(f"integrated lass in fit = {np.trapz(mdtot,t/1e10)*1e10}")

t0plot=1.14622e10 ; nyr=1000000
tt=t0plot+np.arange(nyr)
y=mdeval(tt,g)
plt.plot(tt,y,color="gold")
print(f"integrated lass in finegrid fit = {np.trapz(y,tt/1e10)*1e10}")

fig,ax=plt.subplots(2)
dtmyr=(tt-t0plot)/1e6
mdott=mdeval(tt,g)
ax[0].plot(dtmyr,mdott*1e6)
ax[0].set_xlim(0,1)
mremain=tt*0+1
for i in np.arange(nyr): 
	if i>1: mremain[i]=mremain[i-1]-mdott[i]
ax[1].plot(dtmyr,mremain)
ax[1].set_xlabel(f"Age since {t0plot} (Myr)")
ax[1].set_ylabel("mass")
ax[0].set_ylabel(r"$\dot{{{M}}}$ ($10^{-6} M_\odot$/yr)")
plt.savefig("sun_massloss.pdf")