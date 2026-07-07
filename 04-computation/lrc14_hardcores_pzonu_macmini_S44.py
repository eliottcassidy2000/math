"""
mac-mini-2026-07-07-S44 (HYP-4917) -- THE TWO HARD CORES x PZ-on-U (monad-S3 handoff b).

monad-S3's route census left EXACTLY TWO hard small-parts on the k=8/9 witness-floor
legs: P8 = {9,10,11,12,13} (k=8 clusters) and P9 = {10,11,12,13} (k=9 clusters), with
"huge G_P headroom" (0.447 / 0.525).  The residual regime: spread>108 clusters with
multi-q CRT-blocking.  This script:

 1. EXACT meas(G_P) for both cores (Fraction interval arithmetic) -- pin monad's numerics.
 2. THE PRODUCT-ROUTE ARITHMETIC (exact): with klein-S156's UNCONDITIONAL k=8 Hunter
    floor mu >= 6/49, independence would give G2 >= meas(G_P8)*6/49; how far below
    m_P = 14249/252252 is that, exactly?  What Hunter-floor value closes it at R = 1?
    What R closes it at 6/49?
 3. ADVERSARIAL G2 at the fixed cores (the real object): descend over clusters E
    (0 in E, k points) minimizing G2 = meas(G_P cap {maxgap_E > 1/7}); seeds include
    CRT-blocking spread shapes.  Report restricted-PZ = E[U 1_G]^2/E[U^2 1_G] (rigorous
    lower bound for G2), R = G2/(mGP*mu), and where the minimum sits.
 4. SHARED-SPECTRUM MECHANISM TEST: |G2 - mGP*mu| vs the difference-frequency shared
    mass  SM(E) = sum_{i<j} sum_{m>=1} |gP_hat(m*(e_j-e_i))|^2  (the THM-579-shaped
    prediction: cross-correlation lives on E's difference set hitting G_P's spectrum).
 5. REVERSAL-SYMMETRIZED MOMENTS (owner item, resolved with a one-line derivation):
    symmetrized-PZ mu >= E[U]^2/(E[U^2] + E[U*U∘refl]) <= plain PZ ALWAYS (since
    U >= 0 => cross-moment >= 0) -- the symmetrization does NOT improve the moment
    bound; its value is search-pruning (THM-639-A).  Report rho_rev = E[U U∘refl]/E[U^2]
    per family (= 1 iff palindrome a.e.; an asymmetry invariant).
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
import random as rnd
rnd.seed(44)

MP = F(14249, 252252)
THETA = 1/7

# ---------- 1. exact meas(G_P) ----------
def meas_GP_exact(P):
    """G_P = {x in [0,1): ||p x|| >= 1/14 for all p in P}; exact complement measure."""
    bad = []
    for p in P:
        w = F(1, 14*p)
        for j in range(p+1):
            lo, hi = F(j,p) - w, F(j,p) + w
            bad.append((max(lo, F(0)), min(hi, F(1))))
    bad.sort()
    merged = []
    for lo, hi in bad:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    tot = sum(hi-lo for lo, hi in merged)
    return 1 - tot

P8 = [9,10,11,12,13]
P9 = [10,11,12,13]
m8 = meas_GP_exact(P8); m9 = meas_GP_exact(P9)
print("=== 1. exact G_P measures at the two hard cores ===")
print(f"meas(G_{{9..13}})  = {m8} = {float(m8):.6f}   (monad numeric 0.447)")
print(f"meas(G_{{10..13}}) = {m9} = {float(m9):.6f}   (monad numeric 0.525)")

print("\n=== 2. the product-route arithmetic (exact) ===")
hunter = F(6,49)
prod8 = m8 * hunter
print(f"k=8: meas(G_P8) * 6/49 = {prod8} = {float(prod8):.6f}  vs m_P = {float(MP):.6f}")
print(f"     shortfall = {float(MP - prod8):+.6f}  ({float((MP-prod8)/MP)*100:.2f}% of m_P)")
mu_needed = MP / m8
print(f"     mu-floor needed at R=1: {mu_needed} = {float(mu_needed):.6f} "
      f"(Hunter 6/49 = {float(hunter):.6f}; gap = {float(mu_needed-hunter):+.6f})")
R_needed = MP / prod8
print(f"     R needed at mu = 6/49: {float(R_needed):.4f}  (klein-S155 target was R >= 0.75)")
print(f"k=9: meas(G_P9) = {float(m9):.4f}; bare k=9 mu-floor = 0 (MISTAKE-122) => "
      f"mu-floor needed at R=1: {float(MP/m9):.6f}")

# ---------- numeric G2 machinery ----------
GRID = 100_000
xs = (np.arange(GRID)+0.5)/GRID
def GP_mask(P):
    ok = np.ones(GRID, dtype=bool)
    for p in P:
        d = np.abs(np.mod(p*xs,1.0)-0.5)
        ok &= (0.5-d) >= 1/14
    return ok
MASK8 = GP_mask(P8); MASK9 = GP_mask(P9)
print(f"\n(grid check: |G_P8| = {MASK8.mean():.6f}, |G_P9| = {MASK9.mean():.6f})")

def gapdata(E):
    ph = np.mod(np.outer(xs, np.array(E,float)),1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    mg = g.max(axis=1)
    U = np.clip(g-THETA,0,None).sum(axis=1)
    return mg, U

def core_stats(E, mask, mgp):
    mg, U = gapdata(E)
    good = mg > THETA
    G2 = float((mask & good).mean())
    mu = float(good.mean())
    UG = U*mask
    EUG = float(UG.mean()); EU2G = float((UG*U).mean())
    pzr = EUG*EUG/EU2G if EU2G>0 else 0.0
    R = G2/(mgp*mu) if mu>0 else float('nan')
    return G2, mu, pzr, R

def crt_blocker(k, spread):
    """random k-point 0-anchored cluster covering all residues mod 2..6, spread ~ target."""
    while True:
        E = [0] + sorted(rnd.sample(range(1, spread), k-2)) + [spread]
        E = sorted(set(E))
        if len(E) != k: continue
        if all(len({e % q for e in E}) == q for q in (5,6)):  # covers 5,6 => 2,3 too
            return E

print("\n=== 3. adversarial G2 at the fixed cores ===")
for (label, P, mask, mgp, k) in (("k=8 core P={9..13}", P8, MASK8, float(m8), 8),
                                  ("k=9 core P={10..13}", P9, MASK9, float(m9), 9)):
    best = (9.9, None)
    seeds = [list(range(k))] + [crt_blocker(k, s) for s in (120, 200, 400)] \
          + [[0]+sorted(rnd.sample(range(1,1000),k-1)) for _ in range(2)]
    for E0 in seeds:
        E = list(E0)
        cur = core_stats(E, mask, mgp)[0]
        for it in range(300):
            i = rnd.randrange(1, k)
            cand = E.copy()
            if rnd.random() < 0.6:
                cand[i] = max(1, cand[i] + rnd.choice([-5,-3,-2,-1,1,2,3,5]))
            else:
                cand[i] = rnd.randint(1, 1200)
            cand = sorted(set(cand))
            if len(cand) != k or cand[0] != 0: continue
            cv = core_stats(cand, mask, mgp)[0]
            if cv < cur: E, cur = cand, cv
        if cur < best[0]: best = (cur, E)
    G2, mu, pzr, R = core_stats(best[1], mask, mgp)
    print(f"{label}: min G2 = {G2:.5f} ({G2/float(MP):.2f}x m_P) at E = {best[1]}")
    print(f"    mu = {mu:.4f}  restricted-PZ = {pzr:.5f} ({pzr/float(MP):.2f}x m_P)  R = {R:.4f}")

# ---------- 4. shared-spectrum mechanism ----------
print("\n=== 4. shared-spectrum test at the k=8 core: |G2 - mGP*mu| vs SM(E) ===")
def gp_hat_sq_at(P, freqs):
    """|gP_hat(n)|^2 for the indicator of G_P at given frequencies, via grid FFT-free dot."""
    g = MASK8.astype(float) if P==P8 else MASK9.astype(float)
    out = {}
    for n in freqs:
        c = np.exp(-2j*np.pi*n*xs)
        out[n] = abs(np.dot(g, c)/GRID)**2
    return out

bankE = [list(range(8)),
         [0,1,2,3,4,5,6,120],
         crt_blocker(8,150), crt_blocker(8,400),
         [0,17,34,51,68,85,102,119],       # 17-dilate AP
         [0,1,2,3,4,5,6,7]]
rows=[]
for E in bankE:
    G2, mu, pzr, R = core_stats(E, MASK8, float(m8))
    diffs = sorted({b-a for a in E for b in E if b>a})
    freqs = sorted({m*d for d in diffs for m in (1,2,3)})[:200]
    hat = gp_hat_sq_at(P8, freqs)
    SM = sum(hat[m*d] for d in diffs for m in (1,2,3) if m*d in hat)
    delta = G2 - float(m8)*mu
    rows.append((E, delta, SM, R))
    print(f"E={str(E):44s} Delta={delta:+.5f} SM={SM:.6f} R={R:.4f}")
dl = np.array([abs(r[1]) for r in rows]); sm = np.array([r[2] for r in rows])
if dl.std()>0 and sm.std()>0:
    print(f"corr(|Delta|, SM) = {np.corrcoef(dl, sm)[0,1]:.3f}")

# ---------- 5. reversal-symmetrized moments ----------
print("\n=== 5. reversal-symmetrized PZ (owner item): one-line result + rho_rev table ===")
print("Derivation: V=(U+U∘refl)/2; {V>0} = {U>0} ∪ refl{U>0} has measure <= 2mu;")
print("PZ(V) => mu >= E[U]^2/(E[U^2]+E[U U∘refl]) <= plain PZ since E[U U∘refl] >= 0.")
print("=> symmetrization does NOT improve the moment bound; value = search-pruning.\n")
FAMS = {
 'AP': list(range(1,14)),
 'monad record': [2,4,6,8,10,11,12,13,14,16,18,20,22],
 'EU-min 3adic (non-palin)': [0,30,36,45,50,54,60,63,70,72,81,90,108],
 'random big': sorted(rnd.sample(range(1,2000),13)),
}
for name,E in FAMS.items():
    mg, U = gapdata(E)
    Ur = U[::-1]   # U(1-x) on the symmetric grid
    rho = float((U*Ur).mean()/ (U*U).mean())
    print(f"  rho_rev({name:26s}) = {rho:.4f}" + ("  (palindrome => 1)" if abs(rho-1)<1e-6 else ""))
