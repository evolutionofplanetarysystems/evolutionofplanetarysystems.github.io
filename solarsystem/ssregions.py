import numpy as np
import matplotlib.pyplot as plt
def reverse(x):
        return x[::-1]

import plotsetup
colors=plt.rcParams["axes.prop_cycle"].by_key()['color']


plt.plot([0.1,1e5],[0,1],color="none") ; plt.xscale("log")
plt.xlabel("Orbit radius (au)")
plt.ylabel("Orbital eccentricity")

a=np.logspace(-.5,4.8,1000)
aj=5.2
i=0.

domanytj=False
for tj in [2,2.98]:
	det=aj/a*( (tj-aj/a)/(2*np.cos(i*np.pi/180)))
	e=(1-aj/a*((tj-aj/a)/(2*np.cos(i*np.pi/180)))**2)**.5
	#e[np.isnan(e)]=0
	e[det<0]=1
	e[np.isnan(e)]=0
	q=a*(1-e)
	if domanytj:plt.plot(a,e,label=r"$T_{\rm J}$"+f"={tj:3.1f}")
	if tj>2: ejfc=e
	if tj==2: ehtc=e
if domanytj:plt.legend()
plt.ylim(0,1)
plt.xlim(a.min(),a.max())

pp=np.array([88,225,365,687,4331,10747,30589,59800])/365
ap=pp**(2/3)
ep=np.array([0.206,.007,.017,.094,.049,.052,.047,.010])
#plt.scatter(ap,ep,facecolor="black")
labels=["Me","V",r"E","M","J","S","U","N"]
amin=.45
for aa,ee,ll in zip(ap,ep,labels):
	if aa>amin:
		plt.text(aa,ee,ll,fontsize=10,style="italic",ha="center",va="center")

pd=np.array([1682,90560,557*365,11408*365,74.7*365])/365
ad=pd**(2/3)
ed=np.array([0.116,.244,.441,0.8496,.966])
#plt.scatter(ad,ed,facecolor="none",edgecolor="black")
labels=["Ceres","Pluto","Eris","Sedna","Halley"]
for aa,ee,ll in zip(ad,ed,labels):
	plt.text(aa,ee,ll,fontsize=10,style="italic",ha="center",va="center")

x=np.append(np.append([1.5,aj],reverse(a[(a>1.5)&(a<aj)])),[1.5])
y=np.append(np.append([0,0],reverse(ejfc[(a>1.5)&(a<aj)])),[0])
plt.fill(x,y,color="brown",alpha=.2)
plt.text(2.2,.3,"Asteroids",rotation=90,ha="center",va="center")

plt.fill([amin,1.5,1.5,amin,amin],[0,0,1,1,0],color="red",alpha=.2)
plt.text(.8,.45,"Near-Earth Asteroids",rotation=90,ha="center",va="center")

plpc=200
alpc=plpc**(2/3)

x=np.append( a[a>=alpc],   [a.max(),alpc])
y=np.append( ejfc[a>=alpc], [1,1])
plt.fill(x,y,color="darkorchid",alpha=.2)
plt.text(59,.95,"LPC",rotation=0,ha="center",va="center")


x=np.append(a[(a>1.73)&(a<alpc)],[a[(a>1.73)&(a<alpc)][-1],1.73])
y=np.append(ehtc[(a>1.73)&(a<alpc)], [1,1])
plt.fill(x,y,color="purple",alpha=.2)
plt.text(8.5,.92,"HTC",rotation=0,ha="center",va="center")

x=np.append(a[(a>1.73)&(a<alpc)],[a[(a>1.73)&(a<alpc)][-1],1.73])
y=np.append(ejfc[(a>1.73)&(a<alpc)], [1,1])
plt.fill(x,y,color="blue",alpha=.2)
plt.text(5.5,.6,"JFC",rotation=0,ha="center",va="center")


an=ap[-1] ; en=ep[-1]
x=np.append([aj,an],reverse(a[(a>aj)&(a<an)]))
y=np.append([0,0],  reverse(ejfc[(a>aj)&(a<an)]))
plt.fill(x,y,color="yellow",alpha=.2)
plt.text(14,.3,"Centaurs",rotation=90,ha="center",va="center")

ak=2000
an=165**(2/3)
x=np.append( np.append([an,ak,ak,2000],reverse(a[(a>an)&(a<2000)])), [an])
ejn=.3
y=np.append( np.append([ejn,ejn,0,0],  reverse(ejfc[(a>an)&(a<2000)])), [ejn])
plt.fill(x,y,color="slateblue",alpha=.2)
plt.text(90,.7,"Scattered Disk",rotation=60,ha="center",va="center")

det=an/a*( (tj-an/a)/(2*np.cos(i*np.pi/180)))
tj=3.5
en=(1-an/a*((tj-an/a)/(2*np.cos(i*np.pi/180)))**2)**.5
xn=np.append(a[(a>an)&(a<ak)],[ak,an,an])
yn=np.append(en[(a>an)&(a<ak)],[1,1,en[(a>an)&(a<ak)].min()])
x=np.append( a[(en>ejn)&(a>an)&(a<ak)],[ak,a[(en>ejn)&(a>an)&(a<ak)].min()])
y=np.append(en[(en>ejn)&(a>an)&(a<ak)],[ejn,ejn])
plt.fill(x,y,color="grey",alpha=.2)
plt.text(500,.5,"Detached",rotation=0,ha="center",va="center")


#=np.append( np.append([an,2000],reverse(a[(a>an)&(a<2000)])), [an])
#y=np.append( np.append([0,0],    reverse(ejfc[(a>an)&(a<2000)])), [0])
x=[an,ak,ak,an,an]
y=[0,0,ejn,ejn,0]
plt.fill(x,y,color="aqua",alpha=.2)
plt.text(200,.15,"Kuiper Belt",rotation=0,ha="center",va="center")




plt.fill([2000,20000,20000,2000,2000],[0,0,1,1,0],color="grey",alpha=.2)
plt.text(6e3,.4,"Inner Oort Cloud",rotation=90,ha="center",va="center")
plt.fill([20000,1e5,1e5,20000,20000],[0,0,1,1,0],color="grey",alpha=.4)
plt.text(6e4,.4,"Oort Cloud",rotation=90,ha="center",va="center")

plt.xlim(amin,1e5)
#plt.ylim(0.01,1) ; plt.yscale("log")

plt.tight_layout()

plt.savefig("ssregions.pdf")
