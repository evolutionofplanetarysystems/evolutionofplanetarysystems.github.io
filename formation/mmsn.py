import numpy as np
import matplotlib.pyplot as plt
def reverse(x):
        return x[::-1]

import plotsetup
colors=plt.rcParams["axes.prop_cycle"].by_key()['color']

p=np.array([88,225,365,687,1898,4331,10747,30589,59800,90560])/365
a=p**(2/3)
m=np.array([0.33,4.87,5.97,.642,0.0009,1898,568,86.8,102,.013])*1e27 # g
#plt.scatter(ap,ep,facecolor="black")
labels=["Me","V","E","M","C","J","S","U","N","P"]


# replace Ceres with total ast belt
m[4]=2.39e24 ; labels[4]='AB'

# replace Pluto with KB
m[-1]=m[2]*.2 ; labels[-1]=' KB' ; a[-1]=45

amin=0.3
amax=50


# add in an Oort Cloud
if amax>1000:
	a=np.append(a,1e4)
	m=np.append(m,m[2])
	labels=np.append(labels,'OC')


rmin=a*0 ; rmax=a*0 ; sigma=a*0
for i in range(len(a)):
	if i==0:
		rmin[i]=amin
	else:
		rmin[i]=(a[i-1]*a[i])**.5
	if i==len(a)-1:
		rmax[i]=a[i]*1.5
	else:
		rmax[i]=(a[i+1]*a[i])**.5
edges=np.append(rmin,rmax[-1])
area=np.pi*(rmax**2-rmin**2) * 1.5e13**2

sigma=m/area # g/cm2
#plt.step(a,sigma,where="mid")
plt.stairs(sigma,edges=edges,color="grey")
for ai,si,li in zip(a,sigma,labels):
	if (ai>amin) & (ai<amax):
		plt.text(ai,si,li,ha="center",va="bottom")

aa=np.logspace(-.4,np.log10(amax),1000)
hayashi=1700/aa**1.5
plt.plot(aa,hayashi,color=colors[0],linewidth=2)

# following desch
desch=343/2*(10/aa)**2.168
plt.plot(aa,desch,color=colors[2],linewidth=2)
aaug=a[5:9]
areaaug=area[5:9]
edges=np.append(rmin[5:9],rmax[5:9][-1])
maug=np.array([1750,1411,843,1032])*m[2]
sigmaaug=maug/areaaug
edges=np.append(2.,edges)
aaug=np.append(2.5,aaug)
sigmaaug=np.append(2960,sigmaaug)
#.stairs(sigmaaug,edges=edges,color="black",linewidth=2)
plt.step(aaug,sigmaaug,where='mid',color="black",linewidth=2)

plt.xscale("log") ; plt.yscale("log")
plt.xlim(amin,aa.max())
plt.ylim(sigma[(a>amin)&(a<amax)].min()*.2,hayashi.max())
plt.ylim(sigma[(a>amin)&(a<amax)].min()*.2,desch.max())
plt.xlabel("Distance from Sun (au)")
plt.ylabel(r"Surface density (g/cm$^2$)")
plt.tight_layout()
plt.savefig("mmsn.pdf")
