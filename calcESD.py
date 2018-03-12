import numpy as np
from scipy.stats import t

def esd_critical(alpha, n, i):
    '''
    Compute the critical value for ESD Test
    '''
    df = n - i - 2
    p = 1 - alpha / ( 2 * (n-i+1) )
    t_ppr = t.ppf(p, df)
    ret = (t_ppr * (n-i)) / np.sqrt( (n - i - 1 + t_ppr**2) * (n-i+1) )
    return(ret)

def esd_critical(alpha, n, i):
    df = n - i - 2
    p = 1 - alpha / ( 2 * (n-i+1) )
    t_ppr = t.ppf(p, df)
    ret = (t_ppr * (n-i)) / np.sqrt( (n - i - 1 + t_ppr**2) * (n-i+1) )
    return(ret)

#
def getOutliers(fracs):
    '''
    Calculate outlying fractions with an ESD test until the first
    10 values have been removed (i.e. will not form more than 10 edges)
    fracs is a dict in the form of { window = (fraction,count) }
    '''
    # Start by collecting all fractions
    y = []
    ks = []
    for k,v in fracs.items():
        ks.append(k)
        y.append(v[0])

#    ks = list(fracs.keys())
    n = len(y)
    alpha = 0.000000000001
    tokeep = None
    y2 = np.array(y)
    ks = np.array(ks)
    outl = {} # dict to collect output

    ## Compute test statistic until r=10 values have been
    ## removed from the sample.
    for i in range(1, 10):
        if np.std(y2) == 0:
            break

        ares = abs(y2 - np.mean(y2)) / np.std(y2)
        Ri = max(ares)
        index = np.where(ares == Ri)
        y2 = np.delete(y2,index)

        ## Compute critical value.
        if Ri > esd_critical(alpha,n,i):
            tokeep = i

    # Values to keep
    # Grab indices of relevant values and build new dict
    if tokeep:
        val1 = abs(y - np.mean(y))
        val2 = sorted(val1,reverse=True)[tokeep-1]
        index = np.where(val1 >= val2)
        index = index[0].tolist()
        n = 0

        # Go back to original dict and collect values of interest
        for ind in index:
            outl[ ks[ind] ] = fracs[ks[ind]]
            n += 1

    return outl
