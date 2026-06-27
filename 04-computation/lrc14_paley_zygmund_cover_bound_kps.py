"""
lrc14_paley_zygmund_cover_bound_kps.py  (kind-pasteur-2026-06-27-S31ag)

CRUX 1 attack via SECOND MOMENT. meas(S7)=P(N=0), N=#empty inner sectors.
Paley-Zygmund: P(N>0) >= E[N]^2/E[N^2]  =>  P(N=0) <= 1 - S1^2/(S1+2 S2),
using ONLY the pairwise moments S1=E[N], S2=E[C(N,2)] (E[N^2]=S1+2 S2).

THM-576: cap_10 = 55/91 is PAIRWISE-EXACT. So test: does the second-moment bound
1 - S1^2/(S1+2 S2) close the cap for the BINDING rows k=8,9,10? If it equals cap_k,
CRUX 1's binding row is PROVED by a pairwise (S1,S2) inequality. Find where it works
and where it fails (expect k=10 closest to pairwise; k=8 needs S3,S4).
"""
import sys, itertools, random
from fractions import Fraction as F
from math import comb

def sector_of(p): return int((p % 1) * 7)

def N_distribution(E):
    """q[t] = meas{exactly t inner sectors empty}, t=0..6, exact."""
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b); q = [F(0)] * 7
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        covered = set(sector_of(e * ((x0 + x1) / 2)) for e in E)
        t = 7 - len(covered)
        if 0 <= t <= 6: q[t] += x1 - x0
    return q

def moments(q):
    S1 = sum(t * q[t] for t in range(7))            # E[N]
    S2 = sum(comb(t, 2) * q[t] for t in range(7))   # E[C(N,2)]
    EN2 = sum(t * t * q[t] for t in range(7))       # E[N^2] = S1 + 2 S2
    return S1, S2, EN2

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91)}

def pz_bound(q):
    S1, S2, EN2 = moments(q)
    if EN2 == 0: return F(0), S1, S2
    return 1 - S1 * S1 / EN2, S1, S2

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(91)
    print("Second-moment (Paley-Zygmund) cover bound vs cap_k:")
    print("  P(N=0) <= 1 - S1^2/(S1+2 S2)  [pairwise only];  want <= cap_k")
    for k in (8, 9, 10):
        cap = CAPS[k]
        configs = [("consec", tuple(range(k)))]
        configs.append(("cap-minimizer-ish", tuple(sorted({0,1,5,7,8,9,10,11,12,13} | set(range(k)))) [:k]))
        for _ in range(300):
            cfg = tuple(sorted([0] + random.sample(range(1, 22), k - 1)))
            configs.append(("rand", cfg))
        max_pN0 = F(0); argmax_pN0 = None
        max_pz = F(0); argmax_pz = None
        pz_ge_cap = 0   # how many configs have PZ bound > cap (PZ fails to close)
        for name, E in configs:
            q = N_distribution(E)
            pN0 = q[0]
            pz, S1, S2 = pz_bound(q)
            if pN0 > max_pN0: max_pN0 = pN0; argmax_pN0 = E
            if pz > max_pz: max_pz = pz; argmax_pz = E
            if pz > cap + F(1, 10**9): pz_ge_cap += 1
        print(f"\nk={k}: cap_k = {cap} = {float(cap):.5f}")
        print(f"  max P(N=0) over configs        = {float(max_pN0):.5f}  (<= cap? {max_pN0 <= cap})  at {argmax_pN0}")
        print(f"  max PZ-bound 1-S1^2/E[N^2]      = {float(max_pz):.5f}  (<= cap? {max_pz <= cap})")
        print(f"  configs where PZ-bound > cap    = {pz_ge_cap}/{len(configs)}  "
              f"({'PZ CLOSES the row' if pz_ge_cap==0 else 'PZ too weak (needs higher moments)'})")
        # the consec value explicitly
        qc = N_distribution(tuple(range(k))); pzc, S1c, S2c = pz_bound(qc)
        print(f"  consec_{k}: P(N=0)={float(qc[0]):.5f}  S1={float(S1c):.4f} S2={float(S2c):.4f}  PZ={float(pzc):.5f}")
