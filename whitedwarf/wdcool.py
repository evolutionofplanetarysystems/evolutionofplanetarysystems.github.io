from astropy.io import ascii
import matplotlib.pyplot as plt
from scipy.interpolate import make_smoothing_spline
import numpy as np

def splsmooth(x,y,xnew):
	spl = make_smoothing_spline(x, y, lam=0.1)
	return spl(xnew)
	
import plotsetup

d=ascii.read("wdcool.dat")
lum1=d['col1']
t04=d['col2']
t06=d['col3']
t08=d['col4']
t1= d['col5']
t06he=d['col6']
t06h=d['col7']

d2=ascii.read("wdcool2.dat")
lum2=-d2['col1']
t054=d2['col2']
t055=d2['col3']
t061=d2['col4']
t060a=d2['col5']
t068=d2['col6']
t077=d2['col7']
t087=d2['col8']
t1v2=d2['col9']

t=np.linspace(0,12,1000)

use_winget=False
if use_winget:
	ts=[t04,t06,t08,t1]
	labs=["0.4","0.6","0.8","1"]
	lums=[lum1,lum1,lum1,lum1]
else:
	ts=[t04,t054,t061,t077,t1v2]
	labs=["0.4","0.54","0.61","0.77","1"]
	lums=[lum1,lum2,lum2,lum2,lum2]
	
for tt,lab,lum in zip(ts,labs,lums):
	plt.plot(t,10**splsmooth(tt,lum,t),label=lab)

plt.legend(title="Mass")
plt.xlabel("Age (Gyr)")
plt.ylim(1e-5,.05)
plt.yscale("log")
plt.ylabel("Luminosity / solar")
plt.xlim(0,12)

plt.tight_layout()
plt.savefig("wdcool.pdf")
