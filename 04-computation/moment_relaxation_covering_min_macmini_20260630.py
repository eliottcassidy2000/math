#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S76 -- THE COVERING-MIN LOWER BOUND IS A MOMENT RELAXATION.

Reframe the inner problem M(S) = max_t min_v ||vt|| as a moment / Lasserre relaxation over
measures on the circle, and show its finite LEVELS are the repo's existing threads:

  |L_S(r)| = integral_0^1 prod_{v in S} (1 - g_v(t)) dt,   g_v(t) = 1_{||vt|| < r}   (danger indicator)
           = sum_{A subseteq S} (-1)^|A| I_A,   I_A = meas( intersection_{v in A} danger_v )   (inclusion-exclusion)

Grouped by level k (Bonferroni / moment hierarchy):
  S_k := sum_{|A|=k} I_A = integral C(c(t), k) dt,   c(t) := #{v : ||vt||<r}  (dangers active at t)
  T_m := sum_{k=0}^m (-1)^k S_k   (level-m truncation);  Bonferroni: T_even >= |L_S| >= T_odd.

  level 1  T_1 = 1 - S_1 = 1 - |S|*2r        <- UNION BOUND / Fejer average (HYP-3785 spectral gap; FAILS, <0)
  level 2  adds S_2 = sum pair-overlaps       <- 2nd-moment floor (HYP-3571) = S75 signed correction (HYP-3787)
  level inf (m=|S|)  = exact |L_S|            <- lazy-cut (HYP-3782), the Positivstellensatz

Plus: the EXTREMAL lonely set is a FEW ATOMS (t*=n/Phi6 + symmetric images, measure zero), so the
lonely measure's Toeplitz moment matrix is LOW-RANK PSD -> the covering-min lower bound is a
FLAT-EXTENSION (Curto-Fialkow) certificate; atomic extremum => no finite level is exact.
"""
import numpy as np
from fractions import Fraction as F
from math import comb, gcd

n = 14
Phi6 = n*n - n + 1                      # 183
r_star = F(n, Phi6)                     # 14/183 = covering-min of the construction
CONSTR = list(range(1, n-1)) + [n*(n-1)]   # {1..12, 182}
AP = list(range(1, n))                     # {1..13}, M = 1/14 (non-covering, but the LRC extremal)

def norm_frac(x):
    """||x|| = distance to nearest integer, exact for Fraction x."""
    fr = x - int(x)
    if fr < 0: fr += 1
    return min(fr, 1-fr)

def Mexact(S):
    Sg = sorted(set(S)); cand = set()
    for a in Sg:
        for b in Sg:
            for d in (a-b, a+b):
                if d > 0:
                    for k in range(1, d): cand.add(F(k, d))
    best = F(0); atoms = []
    for t in cand:
        g = min(norm_frac(v*t) for v in Sg)
        if g > best: best = g; atoms = [t]
        elif g == best: atoms.append(t)
    return best, sorted(atoms)

def covers(S, n=n):
    return all(any(v % q == 0 for v in S) for q in range(2, n+1))

print("="*88)
print("PART 1 -- THE EXTREMAL LONELY SET IS A FEW ATOMS (measure zero) => low-rank moment matrix")
print("="*88)
for name, S in [("construction {1..12,182}", CONSTR), ("AP {1..13}", AP)]:
    M, allc = Mexact(S)
    # lonely set at r=M = argmax of G(t)=min_v||vt|| = the atoms achieving M
    atoms = [t for t in allc]
    # dedupe & sort, describe antipodal (iota: t <-> 1-t) structure
    atoms = sorted(set(atoms))
    iota_paired = all((1-t) in set(atoms) for t in atoms)
    denoms = sorted(set(t.denominator for t in atoms))
    print(f"\n{name}: M={M}={float(M):.5f}  covering={covers(S)}")
    print(f"  |lonely atoms at r=M| = {len(atoms)}  (measure zero)")
    print(f"  atoms = {[str(t) for t in atoms]}")
    print(f"  denominators = {denoms}   iota-symmetric (t<->1-t): {iota_paired}")
    # for AP, check the atoms are units mod n (HYP-3571: (Z/14)* in phi/2 antipodal pairs)
    if name.startswith("AP"):
        as_units = sorted(set((t.numerator * (n//gcd(t.numerator, t.denominator))) % n for t in atoms if t.denominator == n))
        units = [u for u in range(1, n) if gcd(u, n) == 1]
        print(f"  atoms with denom n={n}: numerators mod n = {sorted(set(t.numerator % n for t in atoms if t.denominator==n))}; units(Z/{n})* = {units}")

print()
print("="*88)
print("PART 2 -- THE BONFERRONI / MOMENT HIERARCHY: S_k and truncations T_m (exact-on-grid)")
print("="*88)
G = 1 << 22
t = np.arange(G) / G

def danger_count(S, r):
    """c(t) = # dangers active; vectorized on the grid."""
    rf = float(r)
    c = np.zeros(G, dtype=np.int16)
    for v in S:
        d = np.abs(v*t - np.round(v*t))
        c += (d < rf).astype(np.int16)
    return c

def moment_hierarchy(S, r, maxlevel=8):
    c = danger_count(S, r)
    m = len(S)
    S_k = [float(np.mean([comb(int(cc), k) for cc in [0]])) for k in [0]]  # placeholder, recompute properly below
    # S_k = mean over grid of C(c(t), k). Compute via bincount of c.
    counts = np.bincount(c, minlength=m+1).astype(float) / G   # fraction of grid with exactly c dangers
    cvals = np.arange(len(counts))
    S_k = []
    for k in range(0, maxlevel+1):
        Sk = float(np.sum(counts * np.array([comb(int(cv), k) for cv in cvals])))
        S_k.append(Sk)
    LS = counts[0]                      # |L_S| = fraction with 0 dangers
    # Bonferroni truncations T_m = sum_{k<=m} (-1)^k S_k
    T = []
    acc = 0.0
    for k in range(0, maxlevel+1):
        acc += ((-1)**k) * S_k[k]
        T.append(acc)
    return S_k, T, LS, float(np.mean(c)), int(c.max())

# use r slightly below M so |L_S|>0 (at r=M the lonely set is measure zero atoms)
for name, S, M in [("construction", CONSTR, r_star), ("AP {1..13}", AP, F(1, n))]:
    r = M * F(99, 100)                  # 0.99 * M: lonely set has positive measure
    S_k, T, LS, meanc, maxc = moment_hierarchy(S, r, maxlevel=8)
    print(f"\n{name}  r=0.99*M={float(r):.5f}  |S|={len(S)}  mean #dangers={meanc:.3f}  max #dangers={maxc}")
    print(f"  true |L_S(r)| (level-inf, exact-on-grid) = {LS:.6f}")
    print(f"  S_k (level-k moment = sum of k-fold danger overlaps):")
    print(f"    " + "  ".join(f"S_{k}={S_k[k]:.4f}" for k in range(0, 6)))
    print(f"  Bonferroni truncations T_m = sum_{{k<=m}} (-1)^k S_k  (converge to |L_S|):")
    for m in range(0, 8):
        tag = ""
        if m == 1: tag = "  <- UNION BOUND (level 1): 1-|S|*2r"
        if m == 2: tag = "  <- +pair correlations (2nd moment / HYP-3571 / S75)"
        print(f"    T_{m} = {T[m]:+.5f}{tag}")
    print(f"  level-1 value 1-|S|*2r = {1 - len(S)*2*float(r):+.4f}  (USELESS: <0; the finite union bound cannot certify)")

print()
print("="*88)
print("PART 3 -- TRIG-MOMENT (TOEPLITZ) VIEW + FLAT EXTENSION (Curto-Fialkow)")
print("="*88)
print("The lonely measure mu = uniform on the extremal atoms. Its Fourier moments y_k = (1/|A|) sum_a e^{-2pi i k a}")
print("form a Toeplitz matrix M_K=[y_{j-l}] that is PSD (Bochner) with rank = #atoms (Caratheodory-Toeplitz).")

def toeplitz_moment_rank(atoms, Kmax=40):
    A = np.array([float(a) for a in atoms])
    def y(k):
        return np.mean(np.exp(-2j*np.pi*k*A))
    ranks = []
    for K in range(1, Kmax+1):
        Mk = np.array([[y(j-l) for l in range(K+1)] for j in range(K+1)])
        # PSD check (Hermitian): min eigenvalue >= -eps
        eig = np.linalg.eigvalsh((Mk + Mk.conj().T)/2)
        rk = int(np.sum(eig > 1e-9*max(1, eig.max())))
        ranks.append((K, rk, float(eig.min())))
    return ranks

for name, S in [("construction", CONSTR), ("AP {1..13}", AP)]:
    M, atoms = Mexact(S)
    atoms = sorted(set(atoms))
    ranks = toeplitz_moment_rank(atoms, Kmax=30)
    # find where rank stabilizes (flat extension)
    final_rank = ranks[-1][1]
    stab = next((K for K, rk, _ in ranks if rk == len(atoms)), None)
    minmineig = min(e for _,_,e in ranks)
    print(f"\n{name}: #atoms={len(atoms)}  Toeplitz M_K rank -> stabilizes at {final_rank} (=#atoms? {final_rank==len(atoms)})")
    print(f"  rank(M_K) for K=1,5,10,20,30: " + ", ".join(f"{rk}" for K,rk,_ in ranks if K in (1,5,10,20,30)))
    print(f"  min eigenvalue over K (PSD/Bochner, should be >= ~0): {minmineig:.2e}   flat-extension rank stable from K>={stab}")
    # verify localization: mu is off every danger => mean over atoms of 1_{||va||<r}=0 (r=M means ||va||>=M, boundary)
    off = all(norm_frac(v*a) >= M for a in atoms for v in S)
    print(f"  localization: every atom is lonely for every v in S (||va|| >= M): {off}  (=> mu is a lonely measure)")

print()
print("="*88)
print("PART 4 -- WHY NO FINITE LEVEL IS EXACT: the atom needs degree ~ denominator(t*) (uncertainty)")
print("="*88)
print("Resolving a spike at t* of denominator D needs Fourier degree ~ D (uncertainty / Nyquist).")
print("So the moment relaxation must reach degree ~ Phi6(n) to 'see' the atom. Numerology vs lazy-cut #cuts:")
for nn in [8, 10, 12, 14]:
    P = nn*nn - nn + 1
    tstar = F(nn, P)
    print(f"  n={nn}: Phi6={P}  t*=n/Phi6={tstar} denom={tstar.denominator}  (lazy-cut n=12 closed at ~208 cuts; degree scale ~ Phi6)")
print("\n=> finite moment level = klein's Delsarte LP (HYP-3784, near-tight at bounded V, gapped as V grows);")
print("   the gap is HYP-3785's spectral gap = the moment-hierarchy face of S61 'finite-D cannot close Step 3'.")
print("DONE.")
