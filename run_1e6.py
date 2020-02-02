import commands
import time

commands.getoutput("rm -fr runs")
commands.getoutput("mkdir runs")

vals_to_run = []
for i in range(1000):
    vals_to_run.append((int(1000 + i/10.), 4001 + (i*2) % 20))

assert len(set(vals_to_run)) == len(vals_to_run)
print(vals_to_run[:10])


for val in vals_to_run:
    assert val[1] % 2 == 1

    f = open("tmp.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=example
#SBATCH --partition=shared
#SBATCH --time=0-01:50:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile
cd /home/drubin/multiplanelensing
python final_vals_with_uncs_multiplane.py """ + str(val[0]) + """ """ + str(val[1]) + """ | tail""")
    f.close()
    
    time.sleep(0.25)
    commands.getoutput("sbatch tmp.sh")

    
