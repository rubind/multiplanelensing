import numpy as np
from scipy.stats import scoreatpercentile
import sys
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
from multiplane import get_diff_time
import tqdm
import pickle
import glob


def getscore(vals, label = ""):
    score = scoreatpercentile(vals, [15.8655, 50., 84.1345])
    print(label, score)
    return score

def latex_conf(vals, fmt = "%.1f"):
    score = scoreatpercentile(vals, [15.8655, 50., 84.1345])
    errp = score[2] - score[1]
    errm = score[1] - score[0]

    if (fmt % errp) == (fmt % errm):
        latex = ("$" + fmt + " \pm " + fmt + "$") % (score[1], errp)
    else:
        latex = ("$" + fmt + "^{+" + fmt + "}_{-" + fmt + "}$") % (score[1], errp, errm)
    print(latex)
    return latex


def stats_for_one_gal(these_diff_times, key, bins, label = ""):
    these_diff_times = np.array(these_diff_times)

    print(key, "Median ", np.median(these_diff_times))
    print(key, "std()", np.std(these_diff_times))
    print(key, "Median(log10()) ", np.median(np.log10(these_diff_times)))
    print(key, "std(log10())", np.std(np.log10(these_diff_times)))
    print(key, "cred", scoreatpercentile(these_diff_times, [15.8655, 84.1345]))
    print(key, "min", np.min(these_diff_times))

    good_inds = np.where(these_diff_times < 1e9)
    
    NA, bins, NA = plt.hist(np.log10(these_diff_times[good_inds]), bins = bins, label = label + ": log$_{10} (\Delta t)$=" + latex_conf(np.log10(these_diff_times[good_inds]), fmt = "%.2f"))
    if any(these_diff_times < 0):
        bad_count = float(sum(these_diff_times > 1e9))
        xlim = plt.xlim()
        plt.plot(xlim[0], bad_count, "o", label = key + " Distance Larger\nthan GW 170817\n$\Rightarrow$ No Lensing: P=%.1g" % (bad_count/nsamp), color = 'r')
        plt.xlim(xlim)
        
    plt.yscale('log')
    plt.axvline(np.log10(10 - 1.734), color = 'k', zorder = 2)
    ylim = plt.ylim()
    ylim_geo_mean = np.sqrt(ylim[0]*ylim[1])
    plt.text(np.log10(1.734) + 0.025, ylim_geo_mean, "Observed GW-Photon Delay", rotation = 90, size = 9, va = 'center', ha = 'left')
    xticks = plt.xticks()[0]
    xlim = plt.xlim()
    print("xticks", xticks)
    plt.xticks(xticks, ["$10^{%.0f}$" % item for item in xticks])
    plt.xlim(xlim)
    return bins


all_diff_times = {} #pickle.load(open(sys.argv[1], 'rb'))
for fl in glob.glob("runs/*"):
    these_samps = pickle.load(open(fl, 'rb'))
    for key in these_samps:
        try:
            all_diff_times[key] += these_samps[key]
        except:
            all_diff_times[key] = these_samps[key]

    
for key in all_diff_times:
    nsamp = len(all_diff_times[key])

for key in all_diff_times:
    print(key, all_diff_times[key][:5])
    all_diff_times[key] = [-sum(item) for item in all_diff_times[key]]
    print(key, all_diff_times[key][:5])


all_diff_times_concatenate = []
key_time = []

for key in all_diff_times:
    all_diff_times_concatenate += all_diff_times[key]
    key_time.append([np.median(all_diff_times[key]), key])

key_time.sort()
key_time = key_time[::-1]
gal_keys = [item[1] for item in key_time]

print("key_time", key_time)

bins = stats_for_one_gal(all_diff_times_concatenate, "temp", 100)
plt.close()


plt.figure(figsize = (6, 15))
plt.subplot(5,1,5)
stats_for_one_gal(all_diff_times["All"], "All", bins, label = "Multi-Plane Lensing from All")
plt.xlabel("Geometric Time Delay (seconds, log scale)")
plt.ylabel("Histogram of PDF Samples\n(log scale)")
plt.legend(loc = 'upper right')



for i, key in enumerate(gal_keys[1:]):
    plt.subplot(5,1,1+i)
    stats_for_one_gal(all_diff_times[key], key, bins, label = "Only " + key)
    plt.ylabel("Histogram of PDF Samples\n(log scale)")
    plt.legend(loc = 'upper right')
plt.savefig("multiplane.pdf", bbox_inches = 'tight')
