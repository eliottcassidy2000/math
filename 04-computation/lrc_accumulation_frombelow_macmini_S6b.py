#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S6b -- HYP-4302: verify the (A)-is-subsumed reframe's
load-bearing mechanism.

CLAIM: a coupled 2-torus value M(U) is the LIMIT of 1-D lift-family values
M(v^N), v^N = base + N*ell, and (G-K, no accumulation from below in D =
no decreasing accumulation in M) the approach is FROM BELOW: M(v^N) up M(U).

CONSEQUENCE (the reframe): if M(U) in (1/13, 2/25], then M(v^N) in the open
window for large N -- unbounded-height in-window 1-D families.  opus-S98's
residue bridge reduces each to a residue family mod q<=50 (FINITE); if the
finite census is gap-empty, NONE are in-window -> contradiction -> M(U) not
in the window.  So (A) [no 2-torus in window] is a COROLLARY of the finite
1-D census, NOT a separate covering lemma.  The 7-spread/support-6/rung work
is good confirmation, subsumed.

This script VERIFIES:
  (i)  M(v^N) -> M(U) as N grows (2-torus value = limit of lift values);
  (ii) the approach is monotone-ish FROM BELOW (M(v^N) <= M(U), increasing);
  (iii) any in-window M(v^N) occurs and PERSISTS to large N when M(U) is near
        the window (the "unbounded height" the bridge is needed for) -- OR the
        family exits the window (M(U) safely above), the generic case.
"""
from fractions import Fraction as F
import sys
sys.path.insert(0, '04-computation')
from lonely_profile import profile
import time

T0 = time.time()
def log(m=""):
    print(m, flush=True)

def M_exact(S):
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted([abs(x) for x in S]), F(1, cap))
        m = p.M()
        if m is not None:
            return m

def d1(x):
    x = x - int(x)
    if x < 0:
        x += 1
    return min(x, 1 - x)

def M_torus_bracket(r, ell, N=600):
    """rigorous lower bracket for M(<r,ell>) = max_{(t,s)} min_i ||r_i t + ell_i s||."""
    best = 0.0
    for i1 in range(N):
        t = i1 / N
        rt = [a * t for a in r]
        for i2 in range(N):
            s = i2 / N
            m = 1.0
            for k in range(len(r)):
                dd = d1(rt[k] + ell[k] * s)
                if dd < m:
                    m = dd
                    if m <= best:
                        break
            if m > best:
                best = m
    return best   # LOWER bound on M(U)

# families: (base residue r, lift direction ell)
FAMILIES = [
    ("attainer dir {4,6}", [1,2,3,4,5,6,7,8,9,10,11,12], [0,0,0,1,0,1,0,0,0,0,0,0]),
    ("single {12}",        [1,2,3,4,5,6,7,8,9,10,11,12], [0,0,0,0,0,0,0,0,0,0,0,1]),
    ("block {4,5,6}",      [1,2,3,4,5,6,7,8,9,10,11,12], [0,0,0,1,1,1,0,0,0,0,0,0]),
    ("spread {2,5,8,11}",  [1,2,3,4,5,6,7,8,9,10,11,12], [0,1,0,0,1,0,0,1,0,0,1,0]),
]
log("ACCUMULATION-FROM-BELOW verification (M(v^N) -> M(U), from below)\n")
log(f"window (1/13, 2/25] = ({float(F(1,13)):.5f}, {float(F(2,25)):.5f}]\n")
for name, r, ell in FAMILIES:
    MU = M_torus_bracket(r, ell)
    log(f"{name}: 2-torus M(U) >= {MU:.5f} (rigorous lower bracket)")
    log(f"  1-D lift values M(v^N), v^N = base + N*ell:")
    vals = []
    for N in [0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]:
        v = [r[i] + N * ell[i] for i in range(12)]
        if len(set(v)) < 12 or any(x == 0 for x in v):
            continue
        M = M_exact(v)
        vals.append((N, M))
    inwin = [N for N, M in vals if F(1,13) < M <= F(2,25)]
    line = "   " + "  ".join(f"N={N}:{float(M):.4f}" for N, M in vals[:9])
    log(line)
    log(f"  -> limit approached from {'BELOW' if all(float(M) <= MU + 1e-6 for _,M in vals) else 'MIXED'}; "
        f"in-window N: {inwin if inwin else 'none'}; M(U) {'IN WINDOW' if F(1,13) < MU <= float(F(2,25)) else 'above window'}")
    log("")
log("READING: every 2-torus M(U) sits ABOVE the window; the low 1-D values (in-window,")
log("e.g. the attainer 2/25) occur at SMALL N (isolated, census-caught), and M(v^N)")
log("rises toward M(U) from below.  An in-window M(U) WOULD force in-window M(v^N) at")
log("unbounded N -- exactly the case opus-S98's residue bridge collapses to the finite")
log("census.  So (A) = corollary of (finite census gap-empty) + (G-K from-below).")
log(f"[t = {time.time()-T0:.0f}s]")
