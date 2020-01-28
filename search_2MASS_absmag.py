from FileRead import readcol
from time import asctime
from numpy import *
from tqdm import trange
import matplotlib.pyplot as plt
from Spectra import get_zCMB

def get_ang(ra1, dec1, ra2, dec2):
    r1 = array([cos(ra1*pi/180.)*cos(dec1*pi/180.), sin(ra1*pi/180.)*cos(dec1*pi/180.), sin(dec1*pi/180.)], dtype=float64)
    r2 = array([cos(ra2*pi/180.)*cos(dec2*pi/180.), sin(ra2*pi/180.)*cos(dec2*pi/180.), sin(dec2*pi/180.)], dtype=float64)
    
    ang = arccos(dot(r1, r2))*180./pi
    return ang


print(asctime())
[ID, RAdeg, DECdeg, l, b, k_c, h_c, j_c, k_tc, h_tc, j_tc, e_k, e_h, e_j, e_kt, e_ht, e_jt, e_bv, r_iso, r_ext, ba, flgs, NA, ts, v] = readcol("/Users/rubind/Downloads/2mrs_v240/catalog/2mrs_1175_done.dat", 'a,ff,ff,ffffff,ffffff,fffaaaaf')
print(asctime())

zhelio = v/(299792.458)
z = array([get_zCMB(RAdeg[i], DECdeg[i], zhelio[i]) for i in range(len(zhelio))])
#z = zhelio

dist = z*4282.7494
mass = 10.**(-0.4*k_c) * dist**2.


tdelay_ra_dec = []
total_delay = 0.

for i in trange(len(RAdeg)):
    dang = get_ang(197.450375, -23.381486, RAdeg[i], DECdeg[i])
    dang_rad = dang*pi/180.
    impact = dang_rad*dist[i]

    if dist[i] > 1 and dist[i] < 39:
        tdelay = dist[i]*(40 - dist[i]) * (mass[i]/impact)**2.
        total_delay += tdelay
        
        tdelay_ra_dec.append([tdelay, impact, dist[i], mass[i], dang, RAdeg[i], DECdeg[i]])


tdelay_ra_dec.sort()
tdelay_ra_dec = tdelay_ra_dec[::-1]

tdelay_ra_dec = [item + [item[0]/total_delay] for item in tdelay_ra_dec]

print("total_delay", total_delay)
for i in range(20):
    print("tdelay=%.3f impact=%.2f dist=%.1f mass=%.1f dang=%.1f: %.5fd %.5fd frac=%.3f" % tuple(tdelay_ra_dec[i]))
    

plt.hist([item[0] for item in tdelay_ra_dec], bins = 20)
plt.yscale('log')
plt.savefig("all_tdelays.pdf")
plt.close()

