# Adapted from https://github.com/ljmartin/estimating_sem_timeseries/blob/master/sem_utils.py

import numpy as np
from scipy import stats
from scipy import optimize
from scipy.stats import t

import statsmodels.api as sm
from statsmodels.tsa.ar_model import AutoReg


# Original functions from Chodera (adapted from ljmartin)

z95 = 1.959963984540054
def statistical_inefficiency(corr, mintime=3):
    N = corr.size
    C_t = sm.tsa.stattools.acf(corr, fft=True, adjusted=True, nlags=N)
    t_grid = np.arange(N).astype('float')
    g_t = 2.0 * C_t * (1.0 - t_grid / float(N))
    ind = np.where((C_t <= 0) & (t_grid > mintime))[0][0]
    g = 1.0 + g_t[1:ind].sum()
    return max(1.0, g)

def statistical_inefficiency_2(x, mintime=3):
    x = np.asarray(x)
    N = len(x)

    acf = sm.tsa.stattools.acf(
            x,
            fft=True,
            adjusted=True,
            nlags=min(N // 2, 10000)
    )

    t = np.arange(len(acf))
    g_t = 2.0 * acf * (1.0 - t / N)

    negative = np.where((acf <= 0) & (t > mintime))[0]
    if len(negative) == 0:
        ind = len(acf)
    else:
        ind = negative[0]

    g = 1.0 + np.sum(g_t[1:ind])
    return max(1.0, g)

 #and these actually return the values for SEM:
def sem_from_chodera(timeseries):
    autocorrelation_time = statistical_inefficiency(timeseries)
    n = len(timeseries)
    #calculate SEM and return CI:
    sem = np.std(timeseries) / np.sqrt(n/autocorrelation_time)
    return sem

def ci_from_chodera(timeseries):
    sem = sem_from_chodera(timeseries)
    return timeseries.mean()-sem*z95, timeseries.mean()+sem*z95


# --------------------------------------------------------------- 

# Reworked functions



def statistical_inefficiency_reworked(x, mintime=3):
    x = np.asarray(x)
    N = len(x)

    acf = sm.tsa.stattools.acf(
            x,
            fft=True,
            adjusted=True,
            nlags=min(N // 2, 10000)
    )

    t = np.arange(len(acf))
    g_t = 2.0 * acf * (1.0 - t / N)

    negative = np.where((acf <= 0) & (t > mintime))[0]
    if len(negative) == 0:
        ind = len(acf)
    else:
        ind = negative[0]

    g = 1.0 + np.sum(g_t[1:ind])
    return max(1.0, g)

def effective_sample_size_from_g(n, g):
    n_eff = n / g
    return n_eff

def sem_from_chodera_reworked(timeseries):
    autocorrelation_time = statistical_inefficiency_reworked(timeseries)
    n = len(timeseries)
    sem = np.std(timeseries, ddof=1) / np.sqrt(n/autocorrelation_time)
    return sem

def ci_from_chodera_reworked(timeseries):
    sem = sem_from_chodera_reworked(timeseries)
    t95 = t.ppf(0.975, df=effective_sample_size_from_g(len(timeseries), statistical_inefficiency_reworked(timeseries)) - 1)
    return timeseries.mean() - sem * t95, timeseries.mean() + sem * t95




# ---------------------------------------------------------------

# Own functions, statistical mechanics approach 

def autocorrelation_time(series: np.ndarray[float], dt: float) -> float:
    n = len(series)
    series = series - np.mean(series)

    var = np.var(series)
    acf = np.zeros(n)

    for lag in range(n):
        x1 = series[:n - lag]
        x2 = series[lag:]
        acf[lag] = np.mean(x1 * x2) / var

    zeros = np.where(acf < 0)[0]
    cutoff = zeros[0] if len(zeros) > 0 else n

    tau_int = dt * (0.5 + np.sum(acf[1:cutoff]))
    return tau_int

def effective_sample_size(series: np.ndarray[float], dt: float) -> float:
    tau_int = autocorrelation_time(series, dt)
    n = len(series)
    n_eff = n * dt / (2 * tau_int)
    return n_eff

def bootstrap_confidence_interval(series: np.ndarray[float], dt: float, n_bootstrap: int = 1000, confidence_level: float = 0.95) -> tuple[float, float]:
    n_eff = effective_sample_size(series, dt)
    n = len(series)
    indices = np.random.choice(n, size=(n_bootstrap, int(n_eff)), replace=True)
    bootstrap_means = np.array([np.mean(series[idx]) for idx in indices])
    lower_bound = np.percentile(bootstrap_means, (1 - confidence_level) / 2 * 100)
    upper_bound = np.percentile(bootstrap_means, (1 + confidence_level) / 2 * 100)
    return lower_bound, upper_bound

