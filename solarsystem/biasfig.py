import matplotlib.pyplot as plt
import numpy as np

import plotsetup
colors=plt.rcParams["axes.prop_cycle"].by_key()['color']

plt.ion()

npts=100 ; alpha=2

x=10.0**(2*np.arange(npts)/(npts-1))

n=(x/x.max())**(-alpha)

f1=.25 + x*0 # north/south bias

xd=5
f2 = (x>xd)*1 + (x<=xd)*np.exp(-(x-xd)**2/2/(xd/5)**2)

dostep=False

for nn in [n*f1,n*f2]:
	if dostep:
		plt.step(x,nn,where='mid')
	else:
		plt.plot(x,nn)
plt.plot(x,n,color='black')
plt.plot(x,n*f1*f2,'--',color='black')

plt.xscale("log")
plt.yscale("log")
plt.ylim(1,n.max())
plt.xlim(x.min(),x.max())
plt.xlabel("Size of object")
plt.ylabel("Number of objects")

def reverse(x):
	return x[::-1]
	
if ~dostep:
	plt.fill(np.append(x,reverse(x)),np.append(n*f1*f2,n*0+n.min()),color=colors[2],alpha=.1)
	plt.fill(np.append(x,reverse(x)),np.append(n*f1*f2,reverse(n)),color=colors[1],alpha=.1)
#	plt.fill(np.append(x,reverse(x)),np.append(n*f1*f2,reverse(n*f1)),color=colors[1],alpha=.1)
#	plt.fill(np.append(x,reverse(x)),np.append(n*f1,reverse(n)),color=colors[0],alpha=.1)

#plt.text(10,8,r"$\textit{not covered}$",rotation=-33)#,color=colors[1])
plt.text(10,8,"not covered",rotation=-33,fontstyle='italic')#,color=colors[1])
plt.text(1.2,80,"too faint",fontstyle='italic')#,color=colors[1])
plt.text(1.1,2.3,"too faint and not covered",rotation=48,fontsize=10,fontstyle='italic')
plt.text(4,2,"detected",fontstyle='italic')#,color=colors[2])

y=n*f1*f2
plt.scatter(x[::10],y[::10],color="black",zorder=2)

plt.tight_layout()
plt.show()
plt.savefig("biasfig.pdf")