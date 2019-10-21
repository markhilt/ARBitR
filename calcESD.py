import numpy as np
from scipy.stats import t
from scipy.special import logsumexp

def esd_critical(alpha, n, i):
    '''
    Compute the critical value for ESD Test
    '''
    df = n - i - 2
    p = 1 - alpha / ( 2 * (n-i+1) )
    t_ppr = t.ppf(p, df)

    sqrt = np.sqrt( (n - i - 1 + t_ppr**2) * (n-i+1) )
    if sqrt != 0:
        # Some transformations are needed to handle really small sqrt's
        if sqrt < 0.01:
            import ipdb; ipdb.set_trace()
        # np.exp(logsumexp((t_ppr * (n-i))) - logsumexp(sqrt))
        # np.logaddexp(t_ppr * (n-i)) - np.logaddexp(sqrt)
        return (t_ppr * (n-i)) / sqrt
    else:
        return 0


def getOutliers_QC(fracs,windows,r):
    '''Calculate outlying fractions with an ESD test.

    fracs is a numpy array containing fraction values to be searched for outliers.
    windows is a list of the windows in order of appearance in fracs. r is the number of
    values to remove (i.e. number of outliers to consider)
    '''
    # Start by collecting all fractions
    n = len(fracs)
    if r >= n:
        r = n - 2
    alpha = 0.0000000001
    tokeep = None
    #y = fracs.astype("float32")
    y = fracs
    y2 = y
    outl = {} # dict to collect output

    ## Compute test statistic until r values have been
    ## removed from the sample.
    for i in range(1, r):
        if np.std(y2) == 0:
            break

        try:
            ares = abs(y2 - np.mean(y2)) / np.std(y2)
            Ri = max(ares)
            index = np.where(ares == Ri)
            y2 = np.delete(y2,index)
        except:
            import ipdb; ipdb.set_trace()

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

        # Go back to original array and collect values of interest
        for ind in index:
            outl[ windows[ind] ] = fracs[ind]

    return outl

def getOutliers(fracs):
    '''
    Calculate outlying fractions with an ESD test until the first
    10 values have been removed (i.e. will not form more than 10 edges)
    fracs is a dict in the form of { window: fraction }
    '''
    # Start by collecting all fractions
    y = list(fracs.values())
    n = len(y)
    alpha = 0.000000000001
    tokeep = None
    y2 = np.array(y)
    ks = np.array(list(fracs.keys()))
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
