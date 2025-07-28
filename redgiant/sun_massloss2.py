from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

d=ascii.read("0010000M.track.eep")

t=d['star_age']
mdot=-d['star_mdot']

plt.clf()
plt.plot(t,mdot,color="black")
plt.xlim(1.1462e10,1.1463e10)

def mdeval(t):
# t = time since t0=1.14622e10
	dt=89000
	amp=np.array([9.3e-6,8.9e-6,1.5e-6,6e-7,5.2e-7,1.25e-5,1.7e-6,
		2.8e-5,1.3e-6,2.8e-6,
		1.9e-6
		])
	cen=np.array([1.14624100,1.14624115,1.1462413314,1.14623895,1.14626925,1.1462706435,1.14627083,
		1.14628332,1.146281,1.1462834815,
		1.146260018
		])*1e10 - 1.14622e10
	fwhm=np.array([716,862,2000,65000,49000,715,1200,
		900,35000,1760,
		63000
		])
	dmg=1.08*fwhm*amp
	gsum=t*0
	f=(8.*np.log(2))**.5
	for i in range(len(amp)): 
		gi=amp[i]*np.exp(-(t-cen[i])**2/2/(fwhm[i]/f)**2)
		gsum=gsum+gi
	mdlast=.73e-15*np.clip(t-65000,a_min=0,a_max=dt)**2
	mdlast[t>dt]=0
	troll=6000
	delt=t - (dt-troll)
	l=np.where(delt>0)
	rf=t*0+1
	rf[l]=np.exp(-delt[l]**2/.2/troll**2)
	mdlast=mdlast*rf
	mdtot=gsum+mdlast
	return mdtot

t0=1.14622e10
mdtot=mdeval(t-t0)
plt.plot(t,mdtot,color="gray")


print(f"integrated lass in track = {np.trapz(mdot,t/1e10)*1e10}")
print(f"integrated lass in fit = {np.trapz(mdtot,t/1e10)*1e10}")

t0plot=1.14622e10 ; nyr=1000000
tt=t0plot+np.arange(nyr)
y=mdeval(tt-t0)
plt.plot(tt,y,color="gold")
print(f"integrated lass in finegrid fit = {np.trapz(y,tt/1e10)*1e10}")

fig,ax=plt.subplots(2)
dtmyr=(tt-t0plot)/1e6
mdott=mdeval(tt-t0)
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