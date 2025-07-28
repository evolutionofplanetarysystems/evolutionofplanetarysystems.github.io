import matplotlib.pyplot as plt
from runmerc import runmerc
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

yn=input("Run the simulation (y, only if you have MERCURY installed) or just make the plots (return)? ")
if (yn=="y") or (yn=="Y"):
	t,av,ev,iv=runmerc(c_thresh=1e32,runmod=True,object="JUPITER")


thresh=2e31
fig,ax=plt.subplots(3,3)
tmer,amer,emer,imer=runmerc(c_thresh=thresh,runmod=False,object="MERCURY",doplot=True,ax=ax[0,0])
tv,av,ev,iv=runmerc(c_thresh=thresh,runmod=False,object="VENUS",doplot=True,ax=ax[0,0])
te,ae,ee,ie=runmerc(c_thresh=thresh,runmod=False,object="EARTHMOO",doplot=True,ax=ax[0,1])
tm,am,em,im=runmerc(c_thresh=thresh,runmod=False,object="MARS",doplot=True,ax=ax[0,2])
tj,aj,ej,ij=runmerc(c_thresh=thresh,runmod=False,object="JUPITER",doplot=True,ax=ax[1,0])
tsat,asat,esat,isat=runmerc(c_thresh=thresh,runmod=False,object="SATURN",doplot=True,ax=ax[1,1])
tu,au,eu,iu=runmerc(c_thresh=thresh,runmod=False,object="URANUS",doplot=True,ax=ax[1,2])
tnep,anep,enep,inep=runmerc(c_thresh=thresh,runmod=False,object="NEPTUNE",doplot=True,ax=ax[2,0])
tmax=te.max()
plt.tight_layout()
plt.savefig("planet_a.pdf")

#--------------
# Figure of semimajor axis of orbits and radius of star

fig,ax=plt.subplots(figsize=(9,9))
ax.plot(tmer,amer,label="Mercury")
ax.plot(tv,av,label="Venus")
ax.plot(te,ae,label="Earth")
ax.plot(tm,am,label="Mars")
ax.plot(tj,aj,label="Jupiter")
ax.plot(tsat,asat,label="Saturn")
ax.plot(tu,au,label="Uranus")
ax.plot(tnep,anep,label="Neptune")

for t,a,e in zip([tmer,tv,te,tm,tj,tj,tsat,tu,tnep],
	[amer,av,ae,am,aj,aj,asat,au,anep],[emer,ev,ee,em,ej,ej,esat,eu,enep]):
	q=a*(1-e) ; Q=a*(1+e)
	ax.fill(np.append(t,t[-1::-1]), np.append(q,Q[-1::-1]),alpha=.1)
plt.xlim(0,te.max())
ax.ticklabel_format(style='plain', axis='y')
ax.set_yscale("log") 
ax.set_ylabel("a") ; ax.set_xlabel("t")
#ax.ticklabel_format(style='plain', axis='y')
plt.xlabel("Time (years)")
plt.ylabel("a (au)")
plt.ylim(0,)

dstar=ascii.read("/Users/wtreach/Dropbox/EvolutionPlanetary Systems/Figures/track1p0.csv")
tstar=dstar['star_age']-11.4622e9
rstar=10**dstar['log_R']*6.96e10/1.496e13
plt.plot(tstar,rstar,color="grey")
yl=plt.ylim()
ax.fill(np.append(tstar,tstar[-1::-1]),np.append(rstar,rstar*0+yl[0]),
	color="grey",alpha=.2,label="in star")

ax.legend(loc='upper left')
plt.tight_layout()
plt.savefig("allplanet_a.pdf")

#--------------------
# plot of a vs e
plt.figure()
plt.scatter(amer,emer,s=1,label="Mercury")
plt.scatter(av,ev,s=1,label="Venus")
plt.scatter(ae,ee,s=1,label="Earth")
plt.scatter(am,em,s=1,label="Mars")
plt.scatter(aj,ej,s=1,label="Jupiter")
plt.scatter(asat,esat,s=1,label="Saturn")
plt.scatter(au,eu,s=1,label="Uranus")
plt.scatter(anep,enep,s=1,label="Neptune")
plt.legend()
plt.xscale("log")
plt.xlabel("a (au)") ; plt.ylabel("e")
#plt.gca().ticklabel_format(style='plain', axis='x')
plt.savefig("planet_ae.pdf")