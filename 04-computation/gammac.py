"""gammac: finite-R capacity floor gamma_c(R) for AMM 12592 (THM-3029 R1).

RECONSTRUCTED 2026-08-03 (boxeph). The original module lived only in the
deleted session worktree /tmp/math-wt-coinC2/04-computation and was never
committed. The functions below are verbatim from the committed
04-computation/amm12592_finite_R_capacity_threshold.py, with ONE historical
difference recovered from the committed expected output
05-knowledge/results/amm12592_floor_rate_attained_thm3029.out: the default
bisection bracket here is lo=0.5 (the module-level demo table therefore
clamps gamma_c at 1/2 for R = 8, 16, 32 -- exactly what lines 12-20 of that
.out show), whereas THM-3029's referee passes lo=0.02, hi=0.75 explicitly
and so reproduces the unclamped values 0.375, 0.4412, 0.5, ... of the
committed amm12592_finite_R_capacity_threshold.out.

gamma_c(R) = least rate gamma at which the (ARCH) necessary condition holds
at epoch size R (Bernstein capacity criterion of THM-3002 (4) / THM-3027 (*),
evaluated on the extremal shell profile in log space):
    (ARCH) at shell m, extremal profile a_k = min(m-1-k, gamma(m+k)),
    kink L_k: binom(m-1,d) <= sum_k binom(a_k, d-L_k) 2^{a_k-d+L_k} for all d.
Monotone in gamma, so bisection converges to the threshold from above.
GSTAR = lim_R gamma_c(R) = 2 log(phi)/log 5 = log(phi)/log(sqrt 5) (THM-3027).
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

def gamma_c(m, lo=0.5, hi=0.75, iters=44):
    if arch_margin(m, hi) < 0: return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        if arch_margin(m, mid) >= 0: hi = mid
        else: lo = mid
    return hi

# ---------------------------------------------------------------------------
# Module-level demonstration (as in the original lost module): lines 8-25 of
# the committed expected output of THM-3029's referee. Runs at import time.
# ---------------------------------------------------------------------------
print("gamma* (proved asymptotic floor rate) = %.16f" % GSTAR)
print("3/5                                   = 0.6000000000000000")
print()
print("  R      gamma_c(R)      gamma_c - gamma*    is 3/5 admissible?   room below 3/5")
for _R in (8, 16, 32, 64, 128, 256, 512, 1024):  # gamma_c rises to gamma* FROM BELOW
    _gc = gamma_c(_R)
    if _gc is None:
        print(f"{_R:5d}   (no admissible gamma <= 0.75)"); continue
    _ok35 = "YES" if arch_margin(_R, 0.6) >= 0 else "no"
    _room = 0.6 - _gc
    print(f"{_R:5d}   {_gc:.10f}   {_gc-GSTAR:+.10f}      {_ok35:3s}                 {_room:+.10f}")
print()
print("Reading: gamma_c(R) is the NECESSARY-condition floor at finite R.  Any")
print("construction rate must be >= it.  If the empirical construction threshold")
print("gamma_R equals gamma_c(R), the construction is TIGHT against the necessary")
print("condition; if gamma_R is stuck at 3/5 while gamma_c(R) is well below, the")
print("gap is a CONSTRUCTION deficiency, not a capacity one.")
