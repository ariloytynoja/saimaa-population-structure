import sys
import msprime
import numpy as np
import math
import pandas as pd

def simulatePops(npop,cur_size,mig_rate):
    N_a    = 30000          # ancestral size
    N_x    = 1000           # colonization size
    N_c    = cur_size       # current size
    N_p    = int(N_c/npop)  # current subpopulation size

    gen_time       = 10                  # generation time
    T_pop_split    = int(6000/gen_time)  # split time
    T_size_change  = int(10000/gen_time) # colonizing lineage

    gr_p = (math.log(N_p)-math.log(N_x/npop))/T_pop_split

    # subpopulations
    demography = msprime.Demography.island_model([N_p] * npop, 0)
    
    # star-shape migration
    for i in range(1,npop):
        demography.set_migration_rate(0,i,mig_rate)
        demography.set_migration_rate(i,0,mig_rate)
        
    # growth rate
    for i in range(npop):
        demography.add_population_parameters_change(time=0, initial_size=N_p, growth_rate=gr_p, population=i)

    # colonizing population X
    demography.add_population(name="X", initial_size=N_x)
    # split to N subpopulations
    demography.add_population_split(time=T_pop_split, derived=[i for i in range(0, npop)], ancestral="X")
    # size decrease on colonizing X 
    demography.add_population_parameters_change(time=T_size_change, initial_size=N_a, growth_rate=0, population=npop)
    
    return demography

def computeStats(ts,npop,breaks):
    tns = len(ts.samples())
    stp = int(tns/npop)    

    dvt = ts.diversity(windows=breaks).tolist() # diversity total
    dvp = []                                    # diversity subpops

    for i in range(0,tns,stp):
        S = ts.samples()[i:(i+stp)]
        dvp.append(np.ravel(ts.diversity(sample_sets=[S],windows=breaks)).tolist())
    
    return np.vstack((dvt,dvp))


length   = 3750000  # average contig length
ncontigs = 330      # comparable to empirical data

cur_size = 180
smp_size = 120
mut_rate = 2e-8 
rec_rate = 1e-8

breaks   = [x for x in range(0, (length+1), 150000)] # window size (smaller than empirical as no sites are masked)

seed     = None

for npop in [1, 2, 3, 4, 5]:
    for mig_rate in [0.1, 0.05, 0.01, 0.005, 0.001]:

        res = np.empty(shape=(npop+1,ncontigs*25))

        for i in range(0,ncontigs):
       
            demography = simulatePops(npop,cur_size,mig_rate)
            samples = {i: int(smp_size/npop) for i in range(0,npop)}

            ts = msprime.sim_ancestry(samples=samples, demography=demography, sequence_length = length, recombination_rate = rec_rate, random_seed=seed)
            ts = msprime.mutate(ts, rate=mut_rate, model=msprime.InfiniteSites(alphabet=msprime.NUCLEOTIDES))
           
            div = computeStats(ts,npop,breaks) 

            if(i==0):
                res = np.vstack(div)
            else:
                res = np.hstack((res,np.vstack(div)))


        res_df = pd.DataFrame(np.transpose(res), columns = ["tot"]+["pop" + str(x) for x in range(0,npop,1)])
        
        filen  = "runs/"+str(sys.argv[1])+"/simu_"+str(cur_size)+"_"+str(npop)+"_"+str(mig_rate)+".csv"

        res_df.to_csv(filen,index=False)

