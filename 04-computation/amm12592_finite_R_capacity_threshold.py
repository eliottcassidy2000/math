"""Finite-R CAPACITY threshold gamma_c(R): the smallest rate at which the (ARCH)
necessary condition holds at epoch size R.  This is the ceiling any construction
must respect.  Comparing it with the empirical CONSTRUCTION threshold gamma_R
says whether the construction is tight against the necessary condition or is
leaving room on the table.

(ARCH) at shell m, extremal profile a_k = min(m-1-k, gamma(m+k)), kink L_k:
    binom(m-1,d) <= sum_k binom(a_k, d-L_k) 2^{a_k-d+L_k}   for every d.
"""
from math import lgamma, log, exp, sqrt

PHI = (1 + sqrt(5)) / 2
GSTAR = 2 * log(PHI) / log(5)          # 0.5979874356654401497

def lbinom(n, k):
    if k < 0 or k > n or n < 0: return None
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

def lse(v):
    v = [x for x in v if x is not None]
    if not v: return None
    M = max(v)
    return M + log(sum(exp(x - M) for x in v))

def arch_margin(m, gamma):
    """min over d of (log supply - log demand), normalised by m*log2."""
    kappa = (1 - gamma) / (1 + gamma)
    prof = []
    for k in range(m):
        a = min(m - 1 - k, gamma * (m + k))
        L = m * (1 - gamma) - k * (1 + gamma) if k <= kappa * m else 0.0
        if a >= 0 and L >= 0: prof.append((a, L))
    worst = float('inf')
    for d in range(1, m):
        dem = lbinom(m - 1, d)
        if dem is None: continue
        terms = []
        for (a, L) in prof:
            j = d - L
            if j < 0 or j > a: continue
            lb = lbinom(round(a), round(j))
            if lb is None: continue
            terms.append(lb + (a - j) * log(2))
        sup = lse(terms)
        val = ((-1e18) if sup is None else sup) - dem
        if val < worst: worst = val
    return worst / (m * log(2))

def gamma_c(m, lo=0.02, hi=0.75, iters=44):
    if arch_margin(m, hi) < 0: return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        if arch_margin(m, mid) >= 0: hi = mid
        else: lo = mid
    return hi

print("gamma* (proved asymptotic floor rate) = %.16f" % GSTAR)
print("3/5                                   = 0.6000000000000000")
print()
print("  R      gamma_c(R)      gamma_c - gamma*    is 3/5 admissible?   room below 3/5")
for R in (8, 16, 32, 64, 128, 256, 512, 1024):  # gamma_c rises to gamma* FROM BELOW
    gc = gamma_c(R)
    if gc is None:
        print(f"{R:5d}   (no admissible gamma <= 0.75)"); continue
    ok35 = "YES" if arch_margin(R, 0.6) >= 0 else "no"
    room = 0.6 - gc
    print(f"{R:5d}   {gc:.10f}   {gc-GSTAR:+.10f}      {ok35:3s}                 {room:+.10f}")
print()
print("Reading: gamma_c(R) is the NECESSARY-condition floor at finite R.  Any")
print("construction rate must be >= it.  If the empirical construction threshold")
print("gamma_R equals gamma_c(R), the construction is TIGHT against the necessary")
print("condition; if gamma_R is stuck at 3/5 while gamma_c(R) is well below, the")
print("gap is a CONSTRUCTION deficiency, not a capacity one.")
