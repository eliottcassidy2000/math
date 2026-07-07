"""
mac-mini-2026-07-07-S42 -- DIRECT adversarial descent on G2(P,E) = meas(G_P ∩ Good_E)
itself, jointly over (P,E), at every leg k=8..12.

WHY: the hlarge obligation consumes G2 >= m_P. All prior probes bounded G2 from below
via surrogates (union bound needs mu >= thr_k+m_P; monad's condRM = meas(G_P)*(7/6)*
(E[mg|G_P]-1/7) FAILS adversarially at pathology shapes). Nobody has descended on G2
DIRECTLY. If the direct min across joint (P,E) descents stays well above m_P, the
program's residual risk is confirmed to live entirely in (A')-beyond-bands + Part-A
glue, not in a hidden G2 pathology. If some (P,E) approaches m_P, THAT shape is the
real enemy and the fleet needs to see it.

Setup per THM-530: P subset {1..13}, |P| = 13-k; E = cluster co-offsets, 0 in E, |E| = k.
G_P = {x: ||p x|| >= 1/14 for all p in P};  Good_E = {x: maxgap{frac(e x)} > 1/7}.
Includes the THM-530 2/7-pathology shapes as seeds (the known anti-correlated family).
"""
import numpy as np
import random as rnd
from itertools import combinations
rnd.seed(4242)

GRID = 80_000
xs = (np.arange(GRID)+0.5)/GRID
MP = 14249/252252

def G_P_mask(P):
    ok = np.ones(GRID, dtype=bool)
    for p in P:
        d = np.abs(np.mod(p*xs, 1.0) - 0.5)
        ok &= (0.5 - d) >= 1/14   # ||px|| >= 1/14
    return ok

def good_mask(E):
    ph = np.mod(np.outer(xs, np.array(E,float)), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return g.max(axis=1) > 1/7

def G2(P, E):
    return float((G_P_mask(P) & good_mask(E)).mean())

print(f"=== direct G2(P,E) joint adversarial, bar m_P = {MP:.5f} ===")
print("seeds include THM-530's 2/7-pathology shapes (perforated near-AP E, small-heavy P)\n")

overall = []
for k in range(8, 13):
    psz = 13 - k
    best = (9.9, None, None)
    # seed P choices: THM-530 pathology-style (1,2,3-heavy) + spread + random
    seedPs = [tuple(range(1, psz+1))]
    if psz >= 2: seedPs.append(tuple([1,2,3,12,13][:psz]))
    seedPs.append(tuple(sorted(rnd.sample(range(1,14), psz))))
    # seed E choices: perforated near-AP {0,2,3..k} (THM-530 pathology), consecutive, stretched
    seedEs = [tuple([0]+list(range(2, k+1))),
              tuple(range(k)),
              tuple([0]+list(range(2, k))+[2*k])]
    for P0 in seedPs:
        for E0 in seedEs:
            P, E = list(P0), list(E0)
            cur = G2(P, E)
            for step in range(250):
                if rnd.random() < 0.35 and psz > 0:
                    # move P: swap one element
                    cand = P.copy()
                    i = rnd.randrange(psz)
                    cand[i] = rnd.randint(1, 13)
                    if len(set(cand)) != psz: continue
                    cv = G2(cand, E)
                    if cv < cur: P, cur = sorted(cand), cv
                else:
                    # move E: keep 0, perturb another element
                    cand = E.copy()
                    i = rnd.randrange(1, k)
                    if rnd.random() < 0.7:
                        cand[i] = max(1, cand[i] + rnd.choice([-3,-2,-1,1,2,3]))
                    else:
                        cand[i] = rnd.randint(1, 60)
                    cand = sorted(set(cand))
                    if len(cand) != k or cand[0] != 0: continue
                    cv = G2(P, cand)
                    if cv < cur: E, cur = cand, cv
            if cur < best[0]: best = (cur, tuple(P), tuple(E))
    g2v, Pb, Eb = best
    gp = float(G_P_mask(Pb).mean()); mu = float(good_mask(Eb).mean())
    ratio = g2v/MP
    R = g2v/(gp*mu) if gp*mu > 0 else float('nan')
    overall.append((k, g2v, ratio))
    print(f"k={k:2d} (|P|={psz}): min G2 = {g2v:.4f} ({ratio:5.2f}x m_P)  at P={Pb} E={Eb}")
    print(f"        meas(G_P)={gp:.4f}  mu(E)={mu:.4f}  indep product={gp*mu:.4f}  corr ratio R={R:.3f}")

worst = min(overall, key=lambda t: t[1])
print(f"\nWORST direct G2 across all legs: {worst[1]:.4f} at k={worst[0]} = {worst[2]:.2f}x m_P")
print("(numeric probe, grid 80k, 250-step descents x 9 seed pairs per k)")
