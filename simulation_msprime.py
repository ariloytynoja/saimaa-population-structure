import msprime
import numpy as np
import math
import pandas as pd

########################################################################

def simulateOne(length,nind,nseq,mig_rate,debug=False):

    mutation_rate=1e-8
    recombination_rate=1e-8
    gen_time=10
    npop=1

    T_pop_split    = 4000/gen_time
    T_size_change  = 8000/gen_time
    
    N_e    = 10000
    N_a    = 1000
    N_p    = nind

    sample_per_pop = int(nseq/npop)

    gr_p = (math.log(N_p)-math.log(N_a/npop))/T_pop_split
    
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_p, growth_rate=gr_p)]
    
    demographic_events = [      
        msprime.PopulationParametersChange(time=T_pop_split, growth_rate=0, population_id=0),
        msprime.PopulationParametersChange(time=T_size_change, initial_size=N_e, population=0)
    ]
    
    dd = msprime.DemographyDebugger(
            Ne=N_e,
            population_configurations=population_configurations,
            demographic_events=demographic_events)

    if debug:
        dd.print_history() #can comment out when we are happy demography is correct

    ts = msprime.simulate(
        population_configurations = population_configurations,
        demographic_events = demographic_events,
        mutation_rate = mutation_rate,
        length = length,
        recombination_rate = recombination_rate,
        record_migrations = False)
    
    ts = msprime.mutate(
        ts, rate=mutation_rate,
        model=msprime.InfiniteSites(alphabet=msprime.NUCLEOTIDES))

    return ts

########################################################################

def simulateMany(npop,length,nind,nseq,mig_rate,debug=False):

    mutation_rate=1e-8
    recombination_rate=1e-8
    gen_time=10

    T_pop_split    = 4000/gen_time
    T_size_change  = 8000/gen_time
    
    N_e    = 10000
    N_a    = 1000
    N_p    = nind/npop
    
    sample_per_pop = int(nseq/npop)
    gr_p = (math.log(N_p)-math.log(N_a/npop))/T_pop_split
    
    population_configurations = []
    for i in range(npop):
        population_configurations.append(msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_p, growth_rate=gr_p))
    
    ## migration matrix
    migration_matrix = np.zeros((npop,npop))
    migration_matrix[0,1:] = mig_rate
    migration_matrix[1:,0] = mig_rate
    migration_matrix = migration_matrix.tolist()
    
    demographic_events = [
        msprime.PopulationParametersChange(time=T_pop_split, initial_size=N_a, population_id=0),
    ]

    for i in range(1,npop):
        demographic_events.append( msprime.MassMigration(time=T_pop_split, source = i, destination = 0, proportion = 1) )

    demographic_events.append( msprime.PopulationParametersChange(time=T_pop_split, growth_rate=0) )
    demographic_events.append( msprime.MigrationRateChange(time=T_pop_split, rate=0) )
    demographic_events.append( msprime.PopulationParametersChange(time=T_size_change, initial_size=N_e, population=0) )
    
    dd = msprime.DemographyDebugger(
            Ne=N_e,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)

    if debug:
        dd.print_history() 

    ts = msprime.simulate(
        population_configurations = population_configurations,
        migration_matrix = migration_matrix,
        demographic_events = demographic_events,
        mutation_rate = mutation_rate,
        length = length,
        recombination_rate = recombination_rate,
        record_migrations = False)
    
    ts = msprime.mutate(
        ts, rate=mutation_rate,
         model=msprime.InfiniteSites(alphabet=msprime.NUCLEOTIDES))

    return ts

########################################################################

def computeStats(ts,npop):
    tns = len(ts.samples())
    stp = int(tns/npop)    

    npp = [] # number of subpops
    spp = [] # this subpop
    hom = [] # homozygosity
    aff = [] # SAF first
    afl = [] # SAF last
    dvp = [] # diversity this subpop
    nsp = [] # number of segregating sites this subpop
    
    dvt = ts.diversity().tolist()                             # diversity total
    nst = ts.segregating_sites(span_normalise=False).tolist() # number of segregating sites total

    dgf = 0  # genetic divergence central-peripheric 
    dgo = 0  # genetic divergence peripheric-peripheric

    if npop == 2:
        spl=np.array(np.split(np.array(np.arange(0, npop*stp)),npop)).tolist()
        pop_div = ts.divergence(sample_sets=spl)
        dgf = dgo = pop_div.tolist()
        
    if npop > 2:
        spl=np.array(np.split(np.array(np.arange(0, npop*stp)),npop)).tolist()
        pairs = [(i, j) for i in range(npop) for j in range(npop)]
        pop_div = ts.divergence(sample_sets=spl,indexes=pairs).reshape((npop, npop))

        dgf = np.sum(np.triu(pop_div,k=1)[1:,1:])/((npop-1)*(npop-2)/2)
        dgo = np.sum(pop_div[1:,0])/(npop-1)

    for i in range(0,tns,stp):
        S = ts.samples()[i:(i+stp)]
        afs = ts.allele_frequency_spectrum(sample_sets=[S],
                                         polarised=True, span_normalise=False)
        ns = len(afs)-1
        ca = (ns*(ns-1))/2
        sh = 0
        for x in range(2,ns+1):
            a = afs[x]
            ch = (x*(x-1))/2
            sh += a*ch/ca

        npp.append(npop)
        spp.append(int((i+stp)/stp))
        hom.append(sh)
        aff.append(afs[0])
        afl.append(afs[ns])
        dvp.append(ts.diversity(sample_sets=[S]).tolist()[0])
        nsp.append(ts.segregating_sites(sample_sets=[S],span_normalise=False).tolist()[0])

    return ([npp,spp,hom,aff,afl,dvp,nsp],[[npop,dvt,nst,dgf,dgo]])

########################################################################

length = int(100000000)
nind = int(240)
nseq = int(120)
nrep = 100

ar1 = np.array([])
ar1.shape = (8,0)

ar2 = np.array([])
ar2.shape = (8,0)

ar3 = np.array([])
ar3.shape = (8,0)

ar4 = np.array([])
ar4.shape = (8,0)

ar5 = np.array([])
ar5.shape = (8,0)

br1 = np.array([])
br1.shape = (6,0)

br2 = np.array([])
br2.shape = (6,0)

br3 = np.array([])
br3.shape = (6,0)

br4 = np.array([])
br4.shape = (6,0)

br5 = np.array([])
br5.shape = (6,0)


for m in [0, 1e-3, 1e-2, 1e-1]:
    ## 1
    ars = np.array([])
    ars.shape = (7,0)
    brs = np.array([])
    brs.shape = (5,0)
    npop = 1
    for i in range(0,nrep):
        print(m,npop,i, flush= True)
        ts = simulateOne(length,nind,nseq,m)
        (ar,br) = computeStats(ts,npop)
        ars = np.hstack((ars,np.array(ar)))
        brs = np.hstack((brs,np.transpose(np.array(br))))

    mr = np.array([m] * npop*nrep)
    mr.shape = (1,npop*nrep)
    ars = np.vstack((ars,mr))

    mr = np.array([m] * nrep)
    mr.shape = (1,nrep)
    brs = np.vstack((brs,mr))

    ar1 = np.hstack((ar1,ars))
    br1 = np.hstack((br1,brs))

    ## 2
    ars = np.array([])
    ars.shape = (7,0)
    brs = np.array([])
    brs.shape = (5,0)
    npop = 2
    for i in range(0,nrep):
        print(m,npop,i, flush= True)
        ts = simulateMany(npop,length,nind,nseq,m)
        (ar,br) = computeStats(ts,npop)
        ars = np.hstack((ars,np.array(ar)))
        brs = np.hstack((brs,np.transpose(np.array(br))))

    mr = np.array([m] * npop*nrep)
    mr.shape = (1,npop*nrep)
    ars = np.vstack((ars,mr))

    mr = np.array([m] * nrep)
    mr.shape = (1,nrep)
    brs = np.vstack((brs,mr))

    ar2 = np.hstack((ar2,ars))
    br2 = np.hstack((br2,brs))

    ## 3
    ars = np.array([])
    ars.shape = (7,0)
    brs = np.array([])
    brs.shape = (5,0)
    npop = 3
    for i in range(0,nrep):
        print(m,npop,i, flush= True)
        ts = simulateMany(npop,length,nind,nseq,m)
        (ar,br) = computeStats(ts,npop)
        ars = np.hstack((ars,np.array(ar)))
        brs = np.hstack((brs,np.transpose(np.array(br))))

    mr = np.array([m] * npop*nrep)
    mr.shape = (1,npop*nrep)
    ars = np.vstack((ars,mr))

    mr = np.array([m] * nrep)
    mr.shape = (1,nrep)
    brs = np.vstack((brs,mr))

    ar3 = np.hstack((ar3,ars))
    br3 = np.hstack((br3,brs))

    ## 4
    ars = np.array([])
    ars.shape = (7,0)
    brs = np.array([])
    brs.shape = (5,0)
    npop = 4
    for i in range(0,nrep):
        print(m,npop,i, flush= True)
        ts = simulateMany(npop,length,nind,nseq,m)
        (ar,br) = computeStats(ts,npop)
        ars = np.hstack((ars,np.array(ar)))
        brs = np.hstack((brs,np.transpose(np.array(br))))

    mr = np.array([m] * npop*nrep)
    mr.shape = (1,npop*nrep)
    ars = np.vstack((ars,mr))

    mr = np.array([m] * nrep)
    mr.shape = (1,nrep)
    brs = np.vstack((brs,mr))

    ar4 = np.hstack((ar4,ars))
    br4 = np.hstack((br4,brs))

    ## 5
    ars = np.array([])
    ars.shape = (7,0)
    brs = np.array([])
    brs.shape = (5,0)
    npop = 5
    for i in range(0,nrep):
        print(m,npop,i, flush= True)
        ts = simulateMany(npop,length,nind,nseq,m)
        (ar,br) = computeStats(ts,npop)
        ars = np.hstack((ars,np.array(ar)))
        brs = np.hstack((brs,np.transpose(np.array(br))))

    mr = np.array([m] * npop*nrep)
    mr.shape = (1,npop*nrep)
    ars = np.vstack((ars,mr))

    mr = np.array([m] * nrep)
    mr.shape = (1,nrep)
    brs = np.vstack((brs,mr))

    ar5 = np.hstack((ar5,ars))
    br5 = np.hstack((br5,brs))


ars = np.vstack((
         np.transpose(ar1),
         np.transpose(ar2),
         np.transpose(ar3),
         np.transpose(ar4),
         np.transpose(ar5)
    ))

dfA = pd.DataFrame(ars,
  columns=["npp","spp","hom","aff","afl","dvp","nsp","mgr"])

brs = np.vstack((
         np.transpose(br1),
         np.transpose(br2),
         np.transpose(br3),
         np.transpose(br4),
         np.transpose(br5)
    ))

dfB = pd.DataFrame(brs,
  columns=["npp","dvt","nst","dgf","dgo","mgr"])

dfA.to_csv("simulation_resA.csv",index=False)
dfB.to_csv("simulation_resB.csv",index=False)
