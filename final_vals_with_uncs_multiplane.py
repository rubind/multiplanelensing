import numpy as np
from scipy.stats import scoreatpercentile
import sys
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
from multiplane import get_diff_time
import tqdm
import pickle

try:
    from commands import getoutput
except:
    from subprocess import getoutput

def ra_dec_dist_to_xyz(ra, dec, dist):
    return [dist*np.cos(ra*np.pi/180.)*np.cos(dec*np.pi/180.), dist*np.sin(ra*np.pi/180.)*np.cos(dec*np.pi/180.), dist*np.sin(dec*np.pi/180.)]


nsamp = int(sys.argv[1])
Shapiro_steps = int(sys.argv[2])

assert Shapiro_steps % 2 == 1

gwdist = 40.70*10.**(0.2*np.sqrt(0.08**2. + 0.1**2.)*np.random.normal(size = nsamp))
galdists = {"NGC 5084": 16.70*10.**(0.2*0.54*np.random.normal(size = nsamp)),
            "M104": 9.46*10.**(0.2*0.08*np.random.normal(size = nsamp))}

galcoords = {"NGC 5084": [200.070500, -21.827583],
             "M104": [189.997633, -11.623054]}


galmasses = {"NGC 5084": 10**(0.72 + 0.24*np.random.normal(size = nsamp)),
             "M104": 2. * np.exp(0.1*np.random.normal(size = nsamp))}

galxyzs = {}

if 1:
    galdists["NGC 5128"] = 3.76 + 0.05*np.random.normal(size = nsamp)
    galdists["M83"] = 4.79 + 0.1*np.random.normal(size = nsamp)

    galcoords["NGC 5128"] = [201.365063, -43.019112]
    galcoords["M83"] = [204.253958, -29.865417]

    galmasses["NGC 5128"] = 7.2 * np.exp(0.12*np.random.normal(size = nsamp))
    galmasses["M83"] = 0.85 * np.exp(0.06*np.random.normal(size = nsamp))

    # THE HUBBLE FLOW AROUND THE CENTAURUS A/M83 GALAXY COMPLEX
else:
    xyzM83 = [ra_dec_dist_to_xyz(ra = 204.253958, dec = -29.865417, dist = 4.66*10.**(0.2*0.07*np.random.normal())) for i in range(nsamp)]
    xyzCentA = [ra_dec_dist_to_xyz(ra = 201.365063, dec = -43.019112, dist = 3.66*10.**(0.2*0.06*np.random.normal())) for i in range(nsamp)]

    t = np.random.random(size = nsamp)
    galxyzs["Group"] = [np.array(xyzM83[i])*t[i] + np.array(xyzCentA[i])*(1 - t[i]) for i in range(nsamp)]
    galdists["Group"] = [np.sqrt(np.dot(galxyzs["Group"][i], galxyzs["Group"][i])) for i in range(nsamp)]
    galmasses["Group"] = 10. * np.exp(0.25*np.random.normal(size = nsamp))

for key in galmasses:
    galmasses[key] *= galdists[key]/np.median(galdists[key])
    
for key in galcoords:
    galxyzs[key] = [ra_dec_dist_to_xyz(ra = galcoords[key][0], dec = galcoords[key][1], dist = galdists[key][i]) for i in range(nsamp)]



    

all_diff_times = {"All": []}
for key in galmasses:
    all_diff_times[key] = []

for i in tqdm.tqdm(range(nsamp)):
    masses_xyz = []
    masses = []

    for key in galmasses:
        if galdists[key][i] < gwdist[i]:
            masses.append(galmasses[key][i])
            masses_xyz.append(galxyzs[key][i])
            all_diff_times[key].append(get_diff_time(masses_xyz = [galxyzs[key][i]],
                                                     masses = [galmasses[key][i]], source_xyz_init = ra_dec_dist_to_xyz(ra = 197.450375, dec = -23.381486, dist = gwdist[i])))

        else:
            all_diff_times[key].append((-1e10, -1e10))

    diff_time  = get_diff_time(masses_xyz = masses_xyz, masses = masses, source_xyz_init = ra_dec_dist_to_xyz(ra = 197.450375, dec = -23.381486, dist = gwdist[i]), Shapiro_steps = Shapiro_steps)
    print(diff_time)
    all_diff_times["All"].append(diff_time)


getoutput("mkdir runs")
pickle.dump(all_diff_times, open("runs/all_diff_times_%i_steps=%i.pickle" % (nsamp, Shapiro_steps), 'wb'))



"""
elif gal == "5084":
    halovardist = 10**(12.72 + 0.24*np.random.normal(size = nsamp)) 
    dradians = get_ang(197.450375, -23.381486, 200.070500, -21.827583)
    galdist = 16.70*10.**(0.2*0.54*np.random.normal(size = nsamp))
elif gal == "M83":
    log10stellarfixeddist = 11.037 + np.random.normal(size = nsamp)*0.2
    dradians = get_ang(197.450375, -23.381486, 204.253958, -29.865417)
    galdist = 4.66*10.**(0.2*0.07*np.random.normal(size = nsamp))
    galdist0 = 6.960
elif gal == "M104":
    halovardist = 2.e12 * np.exp(0.1*np.random.normal(size = nsamp))
    dradians = get_ang(197.450375, -23.381486, 189.997633, -11.623054)
    galdist = 9.46*10.**(0.2*0.08*np.random.normal(size = nsamp))
elif gal == "CentA":
    halovardist = 1e12
    dradians = get_ang(197.450375, -23.381486, 201.365063, -43.019112)
    galdist = 3.5#*10.**(0.2*0.1*np.random.normal(size = nsamp))

    
print("dradians", dradians)

getscore(gwdist, "GW distance")
latex_conf(gwdist)
getscore(galdist, "Gal distance")
latex_conf(galdist)

try:
    halofixeddist = stellar_to_halo(log10stellarfixeddist)*10.**(np.random.normal(size = nsamp)*0.16)
    getscore(halofixeddist, "Halo, fixed distance")
    getscore(np.log10(halofixeddist), gal + " Log10 Halo, fixed distance")
    latex_conf(np.log10(halofixeddist))

    halovardist = halofixeddist*(galdist/galdist0)**2.
except:
    pass

getscore(halovardist, "Halo, variable distance")

distcombo = galdist*(gwdist - galdist)/gwdist
getscore(distcombo, "Distance Combo")

impact = dradians*galdist
getscore(impact, "Impact (Mpc)")
latex_conf(impact)


timedelay = 18.8682655*(halovardist/1e12)**2. * (distcombo/10.) * impact ** (-2.) # Seconds
getscore(timedelay, gal + " Time delay (s)")
latex_conf(timedelay)
vals = getscore(np.log10(timedelay), "Log10 Time delay (s)")
latex_conf(np.log10(timedelay))


inds = np.where(timedelay > 0)
plt.hist(np.log10(timedelay[inds]), bins = 20)
plt.yscale('log')
plt.title("Time Delay PDF (> 0)")
plt.savefig("time_delay_pdf.pdf", bbox_inches = 'tight')
plt.close()

fracover = sum(timedelay > 1.734 + np.random.normal(size = nsamp)*0.054)/float(nsamp)
print("Fraction over:", fracover)

fracunder = sum(timedelay < 1.734 + np.random.normal(size = nsamp)*0.054)/float(nsamp)
print("Fraction under:", fracunder)


#print("NaN", sum(np.isnan(timedelay)))
#print("Inf", sum(np.isinf(timedelay)))
print("Neg", sum(timedelay < 0))
print("bad distance", sum(galdist > gwdist))

"""
