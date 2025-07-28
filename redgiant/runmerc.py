# Example usage:
# t,a,e,i=runmerc(c_thresh=1e32,runmod=False)
def fixcol(col):
	import numpy as np
	if isinstance(col[0],float):
		res=np.array(col)
	else:
		res=np.zeros(len(col),float)
		for i,c in enumerate(col):
			if isinstance(c,str):
				if c[0]=="*":
					res[i]=np.nan
				else:
					res[i]=float(c)
			else:
				res[i]=c
	return res

def runmerc(runmod=True,object="JUPITER",c_thresh=1e31,pdf=False,doplot=False,ax=None,verbose=True,close=False):
	from os import system
	import numpy as np
	import matplotlib.pyplot as plt
	from astropy.io import ascii
	import time
	import subprocess
	
	if runmod:
		# clean up prior run and run the integrator
		# then convert vectors to elements
		t0=time.time()
		system("rm *tmp *dmp *out")
		system("../mercury6")
		t1=time.time()
		dt=t1-t0
		print(f"=====MERCURY: elapsed time {dt} sec, equal to {dt/60} min or {dt/3600} hr")
		system("rm *aei")
		system("../element6")
		system("rm *.clo")
		system("../close6")
		system("cat info.out")	
		t2=time.time()
		dt2=t2-t0
		print(f"=====RUNMOD: elapsed time {dt} sec, equal to {dt/60} min or {dt/3600} hr")
	
	d=ascii.read(object+".aei",data_start=2)
	time=fixcol(d['col1'])
	long=fixcol(d['col2'])
	M=fixcol(d['col3'])
	a=fixcol(d['col4'])
	e=fixcol(d['col5'])
	incl=fixcol(d['col6'])
	peri=fixcol(d['col7'])
	node=fixcol(d['col8'])
	mass=fixcol(d['col9'])
	q=a*(1-e)
	Q=a*(1+e)
	if doplot:
		ax.plot(time,q,label='perihelion')
		ax.plot(time,a,label='semi-major axis')
		ax.plot(time,Q,label='aphelion')
		ax.set_xlabel("Time (years)")
		ax.set_ylabel("au")
		ax.set_title(object)
		#ax.legend()
	
	if close:	
		close_file=object+".clo"
		command=["cat",close_file]
		res=subprocess.run(command,capture_output=True,text=True)
		lines=len(res.stdout.split('\n'))
		if lines>5:
			dc=ascii.read(close_file,data_start=2)
			ct=fixcol(dc['col1'])
			cp=np.array(dc['col2'])
			cdmin=fixcol(dc['col3'])
			from astropy.io import ascii
			dplanets=ascii.read("planets.ecsv")
			c_score=ct*0
			for i,ccp in enumerate(cp): 
				dd=dplanets[dplanets['planet']==ccp]
				if len(dd)==0: 
					mass=0
				else:
					mass=dd['mass'][0]
				if isinstance(cdmin[i],float):
					c_score[i]=mass/float(cdmin[i])**2
			print("CLOSE ENCOUTNERS")
			print(" time         planet    dmin       score ")
			for i,c_s in enumerate(c_score):
				if c_s>c_thresh:
					print(ct[i],cp[i],cdmin[i],c_s)
					plt.axvline(ct[i],color="black")
					ap=dplanets[dplanets['planet']==cp[i]]['a']
					plt.text(ct[i],ap,cp[i][0])
		else:
			print("no close encounters for "+object)
	if doplot:
		plt.ion();plt.show()
		if pdf:
			plt.savefig("runmod.pdf")
			system("open runmod.pdf")
	return(time,a,e,incl)