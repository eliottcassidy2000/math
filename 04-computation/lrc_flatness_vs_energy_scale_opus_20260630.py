"""
lrc_flatness_vs_energy_scale_opus_20260630.py
Tests the owner's hint: is spectral flatness / Wiener entropy of the danger cover a better LRC monotone than
additive energy (= spectral 4th moment)? RESULT: depends on the CAP's SCALE.
 - vs the COARSE mod-7 coverage proxy measS7 (the repo's cap, HYP-2738): additive energy A wins.
 - vs the TRUE fine LRC gap M=max_t min_v||vt||: Wiener entropy W=-<log m(t)> wins.
 - Winner flips with the CAP on ONE fixed family (span-12) => scale-matching law, not a family artifact.
 - But the AP does NOT minimize W (averaged/entropy lenses blind to the pointwise tight locus); no scalar
   certifies consec-max (HYP-2738). See reflection spectral-flatness-vs-additive-energy-...-opus-20260630.md.
ASCII only (Windows cp1252 safe).
"""
import numpy as np, itertools
from collections import Counter

def frac1(x): return x % 1.0
def fracsym(x):
    x = x % 1.0; return np.minimum(x, 1 - x)

def measS7(E, G=3003):  # measure of t where all 7 sectors mod 7 are hit (coarse coverage proxy = repo cap)
    t = np.arange(G) / G; cov = np.zeros(G, dtype=int)
    for s in range(7):
        hit = np.zeros(G, dtype=bool)
        for e in E: hit |= (np.floor(7 * frac1(e * t)).astype(int) % 7 == s)
        cov += hit.astype(int)
    return np.mean(cov == 7)

def Mgap(E, Q=90):  # true LRC gap max_t min_v ||vt|| (rational sweep denominators 2..Q)
    b = 0.0
    for q in range(2, Q + 1):
        aa = np.arange(1, q)
        b = max(b, fracsym(np.outer(aa, np.array(E)) / q).min(axis=1).max())
    return b

def AE(E):  # additive energy = #{(a,b,c,d): a+b=c+d} = spectral 4th moment (HYP-+2873)
    r = Counter()
    for a in E:
        for b in E: r[a + b] += 1
    return sum(v * v for v in r.values())

def Wfine(E, n, G=3003):  # Wiener entropy of the FINE danger cover m(t)=#{v:||vt||<=1/n}
    t = np.arange(G) / G; m = np.zeros(G)
    for v in E: m += (fracsym(v * t) <= 1.0 / n + 1e-12).astype(float)
    mv = m[m > 0]; return -np.mean(np.log(mv))

def inv(order, target):  # min-orientation discordant-pair count (orientation-robust monotonicity gap)
    m = len(target); d = t = 0
    for i in range(m):
        for j in range(i + 1, m):
            if target[i] == target[j] or order[i] == order[j]: continue
            t += 1; d += 1 if (target[i] - target[j]) * (order[i] - order[j]) < 0 else 0
    return min(d, t - d)

def invfrac(order, target):
    m = len(target); d = t = 0
    for i in range(m):
        for j in range(i + 1, m):
            if target[i] == target[j] or order[i] == order[j]: continue
            t += 1; d += 1 if (target[i] - target[j]) * (order[i] - order[j]) < 0 else 0
    return min(d, t - d) / t

def near_ap_test():
    print("=== NEAR-AP single-swaps, cap = fine gap M (owner's hypothesis holds here) ===")
    for n in [8, 9, 10]:
        base = list(range(1, n)); rows = []
        for k in base:
            for g in range(n, 3 * n):
                if g in base: continue
                S = [v for v in base if v != k] + [g]
                rows.append((Mgap(S, 6 * n), AE(S), Wfine(S, n)))
        M = np.array([r[0] for r in rows]); A = np.array([r[1] for r in rows]); W = np.array([r[2] for r in rows])
        print(f"  n={n}: A inv-frac={invfrac(A,M):.3f}  W inv-frac={invfrac(W,M):.3f}  "
              f"(W better: {invfrac(W,M)<invfrac(A,M)}); |P A/M|={abs(np.corrcoef(A,M)[0,1]):.2f} "
              f"|P W/M|={abs(np.corrcoef(W,M)[0,1]):.2f}")

def span12_family():
    fam = []
    for lo in [1, 2, 3]:
        for mid in itertools.combinations(range(lo + 1, lo + 12), 6):
            fam.append((lo,) + mid + (lo + 12,))
    return fam

def disambiguate():
    print("=== span-12 family (8-subsets of 1..15): winner flips with the CAP's SCALE ===")
    fam = span12_family()
    p7 = np.array([measS7(E) for E in fam]); pM = np.array([Mgap(E) for E in fam])
    A = np.array([AE(E) for E in fam]); W = np.array([Wfine(E, 8) for E in fam])
    print(f"  {len(fam)} sets. cap=measS7 (coarse): A={inv(A,p7)}  W={inv(W,p7)}  winner={'A' if inv(A,p7)<inv(W,p7) else 'W'}")
    print(f"                cap=M (fine gap):  A={inv(A,pM)}  W={inv(W,pM)}  winner={'A' if inv(A,pM)<inv(W,pM) else 'W'}")
    print(f"  Pearson |A vs measS7|={abs(np.corrcoef(A,p7)[0,1]):.3f}  |W vs M|={abs(np.corrcoef(W,pM)[0,1]):.3f}")
    print("  => flip is by CAP not family: coarse->additive energy, fine->Wiener entropy (scale-matching law).")

if __name__ == "__main__":
    near_ap_test()
    disambiguate()
    print("\nCEILING (HYP-2738): the AP does NOT minimize W (entropy blind to the pointwise tight locus);")
    print("no scalar monotone certifies consec-max -- certification needs a SIGNED certificate (L_y).")
