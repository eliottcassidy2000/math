#!/usr/bin/env python3
"""
The covering-min as an ADVERSARIAL FACILITY-LOCATION GAME + a POTENTIAL-FUNCTION /
PRICE-OF-ANARCHY lower bound on the lonely measure (inf meas).  mac-mini-2026-07-01-S91, HYP-3822.

GAME (attacker-defender):
  circle R/Z.  runners = FACILITIES at x_i(t)=v_i t mod 1.  observer = the client at 0.
  loneliness(t) = d(0,t) = min_i ||v_i t||.  M(S)=max_t loneliness = covering objective.
  covering-min = min over covering sets S of M(S) (defender minimizes worst loneliness).
  lonely set L_S(r) = {t : min_i ||v_i t|| >= r}.  inf meas = inf_S |L_S(r)| = the FLOOR.

POTENTIAL (Rosenthal / congestion):  coverage count C(t)=#{i : ||v_i t||<r}.
  A   = int C dt          = (n-1)*2r            (1st moment, total coverage; each |D_i|=2r)
  Phi = int C(C-1) dt      = sum_{i!=j}|D_i cap D_j|  (2nd factorial moment = pairwise congestion)
  C_max = max_t C(t) <= n-1.

THEOREM (potential lower bound on the lonely measure):
  |L_S(r)| >= 1 - A + Phi / C_max.
  PROOF: with m_k=|{C=k}|, A=sum k m_k, |danger|=sum_{k>=1} m_k, Phi=sum k(k-1)m_k.
    A-|danger| = sum_{k>=2}(k-1)m_k;  Phi=sum k(k-1)m_k <= C_max sum(k-1)m_k = C_max(A-|danger|).
    => A-|danger| >= Phi/C_max => |danger|<=A-Phi/C_max => |L|=1-|danger| >= 1-A+Phi/C_max.  QED
  So OVERLAP (congestion) forces loneliness: the union bound (|L|>=1-A) is useless when A>1
  (A~2 at the LRC threshold), but the FORCED congestion Phi>C_max(A-1) makes the floor positive.

PoA reading: the DEFENDER wants to spread facilities (minimize Phi).  The COVERING CONSTRAINT
  (integer speeds, all facilities meet at every rational a/q) forces Phi up = a CONGESTION
  externality = the price of the covering constraint.  inf meas >= 1 - A + inf_S Phi / C_max.
"""
import numpy as np

def danger_sets_and_potential(S, r, G=2_000_000):
    """Fine-grid: coverage count C(t), A, Phi, C_max, |L|(=lonely measure), |danger|."""
    n1 = len(S)
    t = (np.arange(G) + 0.5) / G                      # midpoints of G cells
    C = np.zeros(G, dtype=np.int32)
    for v in S:
        x = (v * t) % 1.0
        d = np.minimum(x, 1.0 - x)                    # ||v t||
        C += (d < r).astype(np.int32)
    A_emp   = C.mean()                                # ~ (n1)*2r
    Phi_emp = (C * (C - 1)).mean()                    # int C(C-1)
    Cmax    = int(C.max())
    L       = (C == 0).mean()                         # lonely measure
    danger  = (C >= 1).mean()
    M_emp   = 0.0                                      # max loneliness (approx via grid min-dist)
    # max_t min_i ||v_i t|| :
    mind = np.full(G, 1.0)
    for v in S:
        x = (v * t) % 1.0
        mind = np.minimum(mind, np.minimum(x, 1.0 - x))
    M_emp = mind.max()
    return dict(A=A_emp, Phi=Phi_emp, Cmax=Cmax, L=L, danger=danger, M=M_emp)

def report(name, S, r, G=2_000_000):
    d = danger_sets_and_potential(S, r, G)
    A, Phi, Cmax, L = d['A'], d['Phi'], d['Cmax'], d['L']
    bound = 1 - A + Phi / Cmax
    union_bound = max(0.0, 1 - A)
    print(f"\n=== {name}  |S|={len(S)}  r={r:.6f} ({r*len(S)*2:.4f}=A_theory) ===")
    print(f"  S = {S}")
    print(f"  A (total coverage, 1st moment) = {A:.5f}")
    print(f"  Phi (Rosenthal congestion pot) = {Phi:.5f}")
    print(f"  C_max (peak coverage, <= n-1)  = {Cmax}")
    print(f"  needed Phi > C_max*(A-1)       = {Cmax*(A-1):.5f}   (Phi {'>' if Phi>Cmax*(A-1) else '<='} it)")
    print(f"  M(S)=max loneliness            = {d['M']:.6f}   (covering-min target 14/183={14/183:.6f})")
    print(f"  --- lonely measure |L| ---")
    print(f"  union-bound floor  1 - A       = {union_bound:.5f}   ({'USELESS' if union_bound<=0 else 'ok'})")
    print(f"  POTENTIAL floor 1-A+Phi/Cmax   = {bound:.5f}   ({'POSITIVE' if bound>0 else 'nonpositive'})")
    print(f"  ACTUAL |L| (lonely measure)    = {L:.5f}")
    print(f"  bound/actual ratio             = {bound/L if L>0 else float('nan'):.3f}")
    return d, bound

if __name__ == "__main__":
    n = 14
    r_thr = 1.0 / n                 # LRC threshold radius 1/14
    r_cov = 14.0 / 183.0            # the covering-min value
    print("#"*72)
    print("# FACILITY-LOCATION POTENTIAL / PoA LOWER BOUND ON THE LONELY MEASURE (n=14)")
    print("#"*72)

    S_ap   = list(range(1, 14))                        # {1..13}
    S_c2   = list(range(1, 13)) + [182]                # {1..12,182} (summary construction)
    # a covering set: multiples enforcing the AP-6 / phi(14) structure -- use {1..13} scaled variants
    S_scal = [ (i* (2)) for i in range(1,14) ]         # {2,4,...,26} (dilation, same M up to scale)

    for r, tag in [(r_thr, "r = 1/14 (LRC threshold)"), (r_cov, "r = 14/183 (covering-min)")]:
        print("\n" + "="*72 + f"\n{tag}\n" + "="*72)
        report("AP {1..13}", S_ap, r)
        report("construction {1..12,182}", S_c2, r)

    # The t=0 forced congestion (the guaranteed pile): Phi_0 >= 4r * sum_{i<j} 1/max(v_i,v_j)
    print("\n" + "#"*72)
    print("# The GUARANTEED (t=0) congestion Phi_0 = 4r * sum_{i<j} 1/max(v_i,v_j)")
    print("#"*72)
    for name, S in [("{1..13}", S_ap), ("{1..12,182}", S_c2)]:
        s = sum(1.0/max(S[i],S[j]) for i in range(len(S)) for j in range(i+1,len(S)))
        for r in (r_thr, r_cov):
            print(f"  {name:14s} r={r:.5f}:  sum 1/max = {s:.4f},  Phi_0 >= {4*r*s:.4f}")
