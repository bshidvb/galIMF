# This code is developed by Li jiadong and Yan Zhiqiang
# Please cite Yan et al. (2024, section 3.1)
# https://ui.adsabs.harvard.edu/abs/2024ApJ...969...95Y/abstract

from matplotlib import rcParams
import matplotlib.pyplot as plt
rcParams["savefig.dpi"] = 100
rcParams["figure.dpi"] = 100
rcParams["font.size"] = 13

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update('jax_platform_name', 'cpu')
from jax import random
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS

### functions:

def Nfrac_lo2hi(N, alpha1=1.3, alpha2=2.3, m_max=150):
    C1 = ((1-alpha1)*N) / (0.5**(1-alpha1) - 0.08**(1-alpha1))
    C2 = C1*(0.5**(alpha2-alpha1))
    return int((C2/(1-alpha2))*(m_max**(1-alpha2)-0.5**(1-alpha2)))

def inv_cumulative_IMF(M, alpha=1.3, m_min=0.1, m_max=0.5):
    '''
    helper function. Just the analytic integral.
    '''
    inside = M*(m_max**(1 - alpha) - m_min**(1 - alpha)) + m_min**(1 - alpha)
    return inside**(1/(1-alpha))

def draw_from_IMF(N, alpha=1.3, m_min=0.08, m_max=150.0):
    '''
    generate N masses from a power-law IMF
    Use inverse-transform sampling, which is much faster than e.g. MC
    '''
    unifs = np.random.uniform(size = N)
    IMF = inv_cumulative_IMF(M = unifs, alpha=alpha, m_min=m_min, m_max=m_max)
    return IMF

def analytic_IMF(M, alpha = 2.35, m_min = 3, m_max = 15):
    '''
    returns the (normalized) IMF, which we can compare to a histogram of masses.
    '''
    c = (1/(1 - alpha) * (m_max**(1 - alpha) - m_min**(1 - alpha)))**(-1)
    return c * M**(-alpha)

def custom_imf_mass_unnormalized(mass, alpha_1, alpha_2, alpha_3):  # there is no time dependence for Kroupa IMF
    if mass < 0.08:
        return 0
    elif mass < 0.5:
        return 0.5**(-alpha_2)/0.5**(-alpha_1)*mass**(-alpha_1+1)
    elif mass < 1:
        return mass**(-alpha_2+1)
    elif mass < 150:
        return mass**(-alpha_3+1)
    else:
        return 0

def mass_function_single_pw(mass, alpha):  # there is no time dependence for Kroupa IMF
    if mass < 0.08:
        return 0
    elif mass < 150:
        return mass**(-alpha) * mass
    else:
        return 0

import jax.numpy as jnp
def pwlaw(mass=None):
    m_min, m_max = jnp.min(mass), jnp.max(mass)
    alpha = numpyro.sample("alpha", dist.Uniform(-5., 5.))

    N = len(mass)
    lgC = numpyro.deterministic(
        "lgC", jnp.log10((jnp.power(m_max, 1. - alpha) - jnp.power(m_min, 1. - alpha)) / (1. - alpha))
    )
    with numpyro.plate("star", N):
        mass_true = mass
        numpyro.factor("lklhd", - lgC - alpha * jnp.log10(mass_true))


### Section 1: random sample stars for a given 2-part power-law IMF
print("This code randomly sample from single power-law IMFs.\n"
      "For a 2-part power-law IMF, it samples stars below and above 0.5 Msun separately.\n"
      "Then it rescale and connect the two samples of stars.\n"
      "Here the plot shows an example of sampling results of the canonical IMF.")

# The total number of stars to be sampled, N_tot
N_tot = 10000
N_lo = 5000  # (this number does not matter)
N_hi = Nfrac_lo2hi(N_lo)

# The rescaled final number of sampled stars in 2 different mass ranges:
N_lo, N_hi = int(N_tot*N_lo/(N_lo+N_hi)), int(N_tot*N_hi/(N_lo+N_hi))

# sample these stars:
mass_lo = draw_from_IMF(N=N_lo, alpha=1.3, m_min=0.1, m_max=0.5)
mass_hi = draw_from_IMF(N=N_hi, alpha=2.3, m_min=0.5, m_max=150)
print("Number of Low-mass sampled stars:", len(mass_lo))
print("Number of High-mass sampled stars:", len(mass_hi))
masses = np.concatenate((mass_lo, mass_hi))  # Join the two arrays
print("Number of joint samples:", len(masses))
# print(masses.shape)

# plot a histogram for the sampling results:
plt.hist(np.log10(masses), bins=50, log=True, histtype='step', label="Joint samples")
plt.hist(np.log10(mass_lo), bins=20, log=True, histtype='step', label="Low-mass sample")
plt.hist(np.log10(mass_hi), bins=30, log=True, histtype='step', label="High-mass samples")
plt.xlabel("Log10 Mass")
plt.legend()
plt.savefig('testplot.pdf', dpi=250)
plt.show()


### Section 2: fit a power law to a sample of (randomly sampled) stars
print("The following code samples a 2-part power-law IMF with variable slopes (alpha1 and alpha2)"
      "Then the code estimate the IMF power-law index for stars between any given stellar mass limits MMIN and MMAX")

# Given any slope for an IMF:
alpha_1 = 1.3
alpha_2 = 2.3
# any upper stellar mass limit:
m_U = 100
# the stars sampled from this given IMF are:
N_tot = 1000000
N_lo = 100000
N_hi = Nfrac_lo2hi(N_lo, alpha1=alpha_1, alpha2=alpha_2, m_max=m_U)
N_lo, N_hi = int(N_tot * N_lo / (N_lo + N_hi)), int(N_tot * N_hi / (N_lo + N_hi))
# print("0.5-{}".format(m_U), N_hi)
masses_rand = []
mass_lo = draw_from_IMF(N=N_lo, alpha=alpha_1, m_min=0.08, m_max=0.5)
mass_hi = draw_from_IMF(N=N_hi, alpha=alpha_2, m_min=0.5, m_max=m_U)
masses = np.concatenate((mass_lo, mass_hi))
masses_rand.append(masses)
masses_rand = np.array(masses_rand).reshape(-1, 1)
# print("0.08-{}".format(m_U), len(masses_rand))
lgmi = np.log10(masses_rand.flatten())

# for any mass limits:
MMIN, MMAX = 0.5, 1
# the stars within the mass limits are:
lgmi = lgmi[(lgmi < np.log10(MMAX)) & (lgmi > np.log10(MMIN))]

print("Summary of inputs: alpha_1=", alpha_1, "alpha_2=", alpha_2, "MMIN=", MMIN, "MMAX=", MMAX)

# use MCMC to fit a single power law to the cut sample.
nuts_kernel = NUTS(pwlaw)
mcmc = MCMC(nuts_kernel, num_samples=800, num_warmup=200)

rng_key = random.PRNGKey(0)
mcmc.run(rng_key, mass=10 ** lgmi)
posterior = mcmc.get_samples()
mcmc.print_summary()

