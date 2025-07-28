import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

import plotsetup
plt.clf()
d=ascii.read("0019000M.track.eep")
plt.plot(d['star_age']/1e9,1.6*1.9/d['star_mass'],label='orbit radius',linestyle="dashed",color="black")
plt.plot(d['star_age']/1e9,10**d['log_R']*6.96e10/1.496e13,label=r"$R_\star$",color="black")
plt.xlim(1.488,1.4917)
plt.ylabel("distance from center of star (au)")
plt.xlabel("age of star (Gyr)")
plt.legend()
plt.tight_layout()
plt.savefig("pollux_rmax.pdf")

