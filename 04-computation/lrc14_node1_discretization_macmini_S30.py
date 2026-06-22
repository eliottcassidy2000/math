"""
Node 1 (THM-527 Part A, finite-Vmax) — discretization lemma + boundary-core closure
(mac-mini-2026-06-22-S30, HYP-2863).

Findings:
1. DISCRETIZATION LEMMA (elementary, rigorous): rho_K = #{good j}/Vmax >= rho* - arcCount/Vmax.
   (Each of the arcCount arcs of GOOD loses <=1 of the Vmax equally-spaced samples to rounding.)
2. BOUNDARY CORE {t,2t,..,12t,V} CLOSED: rho_K depends only on V/t; rho_K>0 for all V/t>12;
   V (observer) > 12t automatic => rho_K>0 always. The s≈0 ruler period collapses the cluster
   teeth to {0,1/2} (maxgap 1/2 > 2/7), a guaranteed good arc of width ~0.065 hit by the j=0
   sample for V/t>7.7. So NO "V/t->inf" needed.
3. Q-UNIFORM: ports to LRC(2q) for all q>=5 (maxgap 1/2 > 2/q; cutoff V/t > 2q-2 = cluster size).
"""
from fractions import Fraction as Fr

def maxgap(pts):
    ph = sorted(set(p % 1 for p in pts))
    if len(ph) <= 1: return Fr(1)
    g = max(b - a for a, b in zip(ph, ph[1:]))
    return max(g, ph[0] + 1 - ph[-1])

def boundary_rhoK(q, t, V):
    """rho_K for the LRC(2q) boundary core {t,..,(2q-2)t, V}, observer V."""
    ncl = 2*q - 2; thr = Fr(2, q)
    coff = [0] + [V - i*t for i in range(1, ncl+1)]
    good = sum(1 for j in range(V) if maxgap([(c*Fr(2*j+1, 2*V)) % 1 for c in coff]) > thr)
    return Fr(good, V)

def discretization_check(cluster, thr, V):
    """Verify rho_K >= rho* - arcCount/V for GOOD={maxgap{frac(e x)}>thr}."""
    E = sorted(set(int(e) for e in cluster)); bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7*abs(e)+1): bp.add(Fr(m, 7*abs(e)))
    bp = sorted(bp); arcs = 0; prev = False; rho = Fr(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        cur = maxgap([(e*((a+b)/2)) % 1 for e in E]) > thr
        if cur: rho += b - a
        if cur and not prev: arcs += 1
        prev = cur
    good = sum(1 for j in range(V) if maxgap([(e*Fr(2*j+1, 2*V)) % 1 for e in E]) > thr)
    return Fr(good, V), rho - Fr(arcs, V), arcs

if __name__ == "__main__":
    print("1. Discretization lemma rho_K >= rho* - arcCount/V (consec_9, thr=2/7):")
    for V in [50, 200, 1000]:
        rk, bound, arcs = discretization_check(list(range(9)), Fr(2, 7), V)
        print(f"   V={V}: rho_K={float(rk):.4f} >= rho*-#arcs/V={float(bound):.4f}: {rk >= bound}")
    print("2. Boundary core rho_K>0 for all V/t>12 (q=7, depends only on V/t):")
    for t in [1, 3, 7]:
        for mult in [13, 15, 20]:
            print(f"   t={t} V/t={mult}: rho_K={float(boundary_rhoK(7, t, mult*t)):.4f}")
    print("3. q-uniform (s≈0 collapse, V/t=2q-1):")
    for q in [5, 7, 9]:
        print(f"   q={q} LRC({2*q}): rho_K={float(boundary_rhoK(q, 2, (2*q-1)*2)):.4f}")
