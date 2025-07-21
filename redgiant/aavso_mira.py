from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

d=ascii.read("aavsodata_mira.csv")
#https://www.aavso.org/data-usage-guidelines#ack
import numpy.ma as nma
 
# remove upper limits
x=np.array(d['Magnitude'])
limit=np.array(np.arange(len(d)),bool)
for i in range(len(x)):
	if x[i][0]=="<":
		limit[i]=True
	else:
		limit[i]=False
dd=d[limit==False]

# select band
dd=dd[dd['Band']=='Vis.']

days=dd['JD']
day0='2015-01-01'
days_since=days-Time(day0).jd
dates=Time(dd['JD'],format='jd')
mags=np.array(dd['Magnitude'],float)
fluxs=3636*10**(-0.4*mags)
fluxrel=fluxs/0.4

plt.figure(figsize=(7.5,3))

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

xd=np.linspace(-1000,2000,100000)
phase=275
amp=300
period=325
x=(xd-phase)*2*np.pi/period
model=amp*10**(np.cos(x)-1)

plt.scatter(days_since,fluxs,3,color="black")
plt.xlim(0,days_since.max())
plt.ylim(0.5,599)
#plt.yscale("log")
plt.xlabel(f"days since {day0}")
plt.ylabel("Brightness / Minimum ")
plt.tight_layout()
#plt.plot(xd,model)
plt.savefig("aavso_mira.pdf")

print("date range from ",day0," through ",dates[-1].iso)

#fig,ax=plt.subplots()
#ax.xaxis.axis_date()
#ax.scatter(dates,mags)
#fig.autofmt_xdate()
