from astropy.table import Table
import glob
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_smoothing_spline

'''
from distance 5.8 au
psg_rad_jup = jupiter at 5.2 au from sun
psg_rad_arcturus = jupiter at 5.2 au from arcturus
psg_rad_mira = jupiter at 5.2 au from mira
psg_rad_wd = jupiter at 10 au from WD (T=1e4, R=0.0097 Rsolar)
# Radiance unit: Spectral irradiance [Jy]
'''

import plotsetup
colors=plt.rcParams["axes.prop_cycle"].by_key()['color']

# make a model spectrum of protostar plus Jupiter...
wp=np.logspace(-1,3,1000)
nup=2.998e14/wp
from planck import planck
tp=300 ; wcross=3.5
from solar_sp import solar_sp
fs=solar_sp(wp)
alb=.08+.7*np.exp(-(wp-.7)**4/2/.3**2)
wa=[1.28,1.6,2.7,1.9] ; aa=[.2,.5,.4,.3]
for waa,aaa in zip(wa,aa): alb=alb+aaa*np.exp(-(wp-waa)**2/2/.05**2)
fs=fs*alb
fp=planck(wp,tp)
lcross=np.argmin(np.abs(wp-wcross))
fc=nup*(fs + fs[lcross]/planck(wcross,tp)*fp)
fc=1.5e23/max(fc) * fc
# Halpha for fun
wline=[.65,.48]
fline=np.array([2,.5])*1e23
for w,f in zip(wline,fline): fc=fc+f*np.exp(-(wp-w)**2/2/.01**2)
lself_now=4.6e21 # W/m2 of current jupiter
fc=fc/np.max(fc)*lself_now*20
plt.plot(wp,fc,color="purple",label="Protostar")


files=glob.glob("psg*.txt")
files.sort()
choose=False
if choose:
	for i,f in enumerate(files): print(f"{i}: {f}")
	sel=input("Which PSG file [0]?")
	if sel=="": 
		isel=0
	else:
		isel=int(sel)
	filesel=files[isel]
	factors=[1]
else:
	files=np.array(files)
	filesel=['psg_jup_Main Sequence.txt','psg_jup_Red Giant.txt','psg_jup_White Dwarf.txt',
		'psg_jup_No Star.txt']
	colors=[colors[2],colors[3],colors[0],'black']
	factors=[1.2,1,.5,.4]
fmin=2e18
fmaxs=[]

# replace white dwarf spectrum with Koester
file='koester_T10000_log8.txt'
d=ascii.read(file)
w=d['col1']/1e4
nu=2.998e14/w
fnu=d['col2']*d['col1']/nu
fnu = fnu/fnu.max()
nufnu=nu*fnu
wwd=w ; fwd=nufnu/nufnu.max()*1.2e-3*1.9e33/1e3*(69911e5/10/1.5e13)**2/2/2

for i,file in enumerate(filesel):
	d=ascii.read(file)
	w=d['col1']
	f=d['col2']
	f=1e-26*f*2.998e14/w
# luminosity
	f=f*4*np.pi*(5.6*1.496e13)**2
	l=np.where(f!=0)
	w=w[l] ; f=f[l]
	#l=np.where(f>fmin)
	#w=w[l] ; f=f[l]

	if file=="psg_jup_White Dwarf.txt":
		wsplice=3.8
		worig=w*1.
		w=np.append(wwd[wwd<wsplice],w[w>wsplice])
		f=np.append(fwd[wwd<wsplice],f[worig>wsplice])
	else:
		f=f*factors[i]
	color=colors[i]

	# smooth
	dosmooth=True
	if dosmooth:
		from scipy.ndimage import gaussian_filter
		#a=np.argsort(w)
		#d=w[1:] - w[:-1]
		#l=np.where(d>0)
		#spl = make_smoothing_spline(w[l], f[l], lam=0.01)
		#f=spl(w)
		f=gaussian_filter(f,3)
	plt.plot(w,f,label=file.split('_')[-1].split('.')[0],color=color)
	fmaxs=np.append(fmaxs,np.nanmax(f))
plt.xscale("log")
plt.yscale("log")
ymax=fmaxs.max()*1.5
plt.ylim(1e16,ymax)
plt.xlim(0.15,500)
plt.ylabel(r"Luminosity ($\nu L_\nu$ in W/m$^2$)")
plt.xlabel(r"Wavelength ($\mu$m)")
plt.legend(loc="lower right",fontsize=10)

plt.text(9,ymax,'Thermal Emission')
plt.text(.7,ymax,'Reflected',ha='center',va='bottom')
plt.tight_layout()
plt.savefig("psg_jup.pdf")

