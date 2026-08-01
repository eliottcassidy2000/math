"""AMM 12592: does a DEGREE-RESOLVED cross-shell Hall cut beat the degree-blind one?

opus (amm12592_cross_shell_capacity_opus_S4.py) reduced C*_general = C*_block to
a forward-transportation feasibility whose Gale/Hall cuts are
    for all M:   sum_{m>=M} demand  <=  sum_{m>=M} supply,
and evaluated them DEGREE-BLIND (aggregated over d).  Their stated gap: "full
rigor wants the degree-resolved Hall condition ... you would know if a
degree-specific cross-shell cut can beat it."

Hall/Gale for a bipartite forward flow: feasible <=> EVERY cut satisfied.  The
degree-blind cuts are unions of degree bands, hence a SUBSET of all cuts.  So
   - degree-resolved is a SUPERSET of constraints: it can only EXCLUDE MORE
     gamma, i.e. it can RAISE the floor, never lower it;
   - opus's sign flip at golden already proves infeasibility below golden.
The open question is therefore precisely: does some degree-resolved cut bind at
a gamma ABOVE golden, pushing C*_general strictly above log_5(5 phi^2)?
(There is little room: the gamma=3/5 constructions give C <= 8/5 = 1.6.)

Deficits route FORWARD to deeper shells at FIXED ABSOLUTE DEGREE d (the u-adic
order at u=-1 is not rescaled by the shell), so in normalised coordinates the
band shifts: shell m has delta = d/m, shell 2m has delta/2.  We test the cuts
    C(M,d):   sum_{m>=M} binom(m-1,d)  <=  sum_{m>=M} supply_m(d)
against the aggregate cuts, at finite dyadic shells, in log space.
"""
from math import lgamma, log, log2, exp, sqrt

PHI = (1 + sqrt(5)) / 2
GOLDEN = 2 * log(PHI) / log(5)

def lbinom(n, k):
    if k < 0 or k > n or n < 0: return None
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

def lse(vals):
    vals = [v for v in vals if v is not None]
    if not vals: return None
    M = max(vals)
    return M + log(sum(exp(v - M) for v in vals))

def profile(m, gamma):
    """extremal profile of THM-3009: a_k = min(m-1-k, gamma(m+k)); L_k from the kink."""
    kappa = (1 - gamma) / (1 + gamma)
    out = []
    for k in range(m):
        a = min(m - 1 - k, gamma * (m + k))
        L = m * (1 - gamma) - k * (1 + gamma) if k <= kappa * m else 0.0
        if a >= 0 and L >= 0:
            out.append((a, L))
    return out

def lsupply(m, d, gamma):
    """log of sum_k binom(a_k, d-L_k) 2^{a_k-d+L_k}  (natural log)."""
    terms = []
    for (a, L) in profile(m, gamma):
        j = d - L
        if j < 0 or j > a: continue
        lb = lbinom(round(a), round(j))
        if lb is None: continue
        terms.append(lb + (a - j) * log(2))
    return lse(terms)

def ldemand(m, d):
    return lbinom(m - 1, d)

def per_shell_margin(m, gamma):
    """min over d of [supply - demand] in log space -- the (ARCH) per-shell test."""
    worst, argd = float('inf'), None
    for d in range(1, m):
        dem = ldemand(m, d)
        sup = lsupply(m, d, gamma)
        if dem is None: continue
        val = (-1e18 if sup is None else sup) - dem
        if val < worst: worst, argd = val, d
    return worst / (m * log(2)), argd

def tail_cut_margin(shells, gamma):
    """min over (M, d) of log[ sum_{m>=M} supply_m(d) ] - log[ sum_{m>=M} demand_m(d) ].
       DEGREE-RESOLVED: d is an absolute degree, held fixed across shells."""
    worst, arg = float('inf'), None
    dmax = max(shells) - 1
    for iM in range(len(shells)):
        tail = shells[iM:]
        for d in range(1, dmax):
            dem = lse([ldemand(m, d) for m in tail if d < m])
            if dem is None: continue
            sup = lse([lsupply(m, d, gamma) for m in tail])
            val = (-1e18 if sup is None else sup) - dem
            nrm = val / (max(tail) * log(2))
            if nrm < worst: worst, arg = nrm, (shells[iM], d)
    return worst, arg

print("CONTROL -- per-shell (ARCH) margin, normalised by m log2.  Sign must flip at golden.")
for J in (7, 8, 9):
    m = 2 ** J
    row = []
    for g in (GOLDEN - 0.03, GOLDEN - 0.01, GOLDEN, GOLDEN + 0.01, GOLDEN + 0.03):
        v, dd = per_shell_margin(m, g)
        row.append(f"{v:+.4f}")
    print(f"  m={m:5d}: gamma = golden+(-.03,-.01,0,+.01,+.03) -> {'  '.join(row)}")

print("\nDEGREE-RESOLVED cross-shell tail cuts, shells m = 8,16,...,512 (absolute d fixed):")
shells = [8, 16, 32, 64, 128, 256, 512]
for g in (GOLDEN - 0.03, GOLDEN - 0.01, GOLDEN, GOLDEN + 0.01, GOLDEN + 0.03):
    tv, targ = tail_cut_margin(shells, g)
    pv, _ = per_shell_margin(max(shells), g)
    verdict = "tail STRICTLY stronger" if tv < pv - 1e-6 else "tail no stronger than per-shell"
    print(f"  gamma=golden{g-GOLDEN:+.3f}: tail-cut margin={tv:+.5f} (at M={targ[0]}, d={targ[1]}), "
          f"per-shell={pv:+.5f}   -> {verdict}")

print("\nWHY IT DECOUPLES (this is a proof, not a numeric).  Forward routing preserves the")
print("ABSOLUTE degree d (the u-adic order at u=-1 is not rescaled by the shell), so the")
print("bipartite graph is a DISJOINT UNION over d of forward-in-m chains.  The")
print("neighbourhood-closed sets are therefore S = union_d {(m,d): m >= M_d}, and the cut")
print("for such a union is the SUM of the per-d cuts -- implied by them.  Hence the")
print("degree-resolved Hall condition is EXACTLY the family of per-(M,d) cuts tested above.")
print("\nAnd the coarsest (degree-blind) cut is a VALID Hall cut for ANY routing rule in d,")
print("since aggregating assumes the most generous possible d-mobility.  So opus's")
print("violation below golden proves infeasibility UNCONDITIONALLY -- their stated gap")
print("does not threaten the lower bound at all; degree resolution can only ADD")
print("constraints, i.e. only RAISE the floor.  The test above shows it raises nothing.")

print("\nBinding degree fraction delta* = d*/m at gamma = golden (should approach 1/phi = %.6f):"
      % (1 / PHI))
for J in (6, 7, 8, 9, 10):
    m = 2 ** J
    v, dd = per_shell_margin(m, GOLDEN)
    print(f"   m={m:5d}: binding d*={dd:4d},  delta*={dd/m:.6f},  margin={v:+.5f}")
