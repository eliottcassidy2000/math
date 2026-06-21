#!/usr/bin/env python3
"""
lrc14_threadB_rotation_vs_paley_symmetry_macmini.py
mac-mini-2026-06-21 -- THREAD B verdict computation.

CANON LEM-004 (scope flags 1&2) already PROVES: Paley is the circulant H-max
ONLY at n=7,11. At n=13,15,17,19 the ROTATION (consec / carousel) tournament
S={1,...,m} is the circulant H-max; at n=19 Paley is 3rd. The DET(I+S)
tournament-Hadamard ceiling (THM-472) is a DRT/spectral extremality, but the
H-max (Hamiltonian-path count) is a STRICTLY FINER, NON-spectral one (THM-499).

THREAD B asks: is the Doyle-Holt / half-arc-transitive (arc-transitive) lens
the H-max characterization? CANON answer is already NO for H. This script makes
the symmetry side EXPLICIT and pins down which extremality each symmetry class
carries, by computing for ROTATION vs PALEY (and tower-DRT) at n=15,19:
  - is it a DRT? (det(I+S) ceiling / weakly-regular spectral extremal)
  - |Aut|, vertex-orbits, arc-orbits, self-converse, arc-transitive
  - H (the actual maximand)
This isolates: the ARC-TRANSITIVE object (Paley) is the det-ceiling / spectral
extremal; the H-max object (rotation) is NOT arc-transitive and NOT a DRT.

n=19 H is expensive (Held-Karp 2^19*19); we use the canon EXACT values from
LEM-004 instead of recomputing, and compute only the (cheap) symmetry data here.
At n=15 we compute H directly (already have tower-DRT H=198335025 from canon).
"""
import sys, time
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

def banner(t):
    print("\n" + "=" * 74 + "\n  " + t + "\n" + "=" * 74)

def circ(n, S):
    S = set(s % n for s in S)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A

def is_tournament(A):
    n = len(A)
    return all(A[i][i] == 0 for i in range(n)) and \
        all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(i + 1, n))

def skewM(A):
    n = len(A)
    return np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=np.int64)

def is_DRT_adj(A):
    n = len(A); S = skewM(A)
    J = np.ones((n, n), dtype=np.int64)
    return np.array_equal(S @ S, J - n * np.eye(n, dtype=np.int64))

def det_I_plus_S(A):
    n = len(A)
    return round(float(np.linalg.det(np.eye(n) + skewM(A))))

def automorphisms(A, cap=None):
    n = len(A)
    outdeg = [sum(A[i]) for i in range(n)]
    by_deg = {}
    for v in range(n):
        by_deg.setdefault(outdeg[v], []).append(v)
    perm = [-1] * n; used = [False] * n
    order = sorted(range(n), key=lambda v: len(by_deg[outdeg[v]]))
    auts = []
    def consistent(i, cand):
        for j in range(n):
            pj = perm[j]
            if pj == -1: continue
            if A[cand][pj] != A[i][j]: return False
            if A[pj][cand] != A[j][i]: return False
        return True
    def bt(k):
        if cap is not None and len(auts) >= cap: return
        if k == n:
            auts.append(perm[:]); return
        v = order[k]
        for cand in by_deg[outdeg[v]]:
            if used[cand]: continue
            if consistent(v, cand):
                perm[v] = cand; used[cand] = True
                bt(k + 1); used[cand] = False; perm[v] = -1
    bt(0)
    return auts

def vertex_orbits(A, auts):
    n = len(A); p = list(range(n))
    def f(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    for g in auts:
        for v in range(n):
            a, b = f(v), f(g[v])
            if a != b: p[a] = b
    return len({f(v) for v in range(n)})

def arc_orbits(A, auts):
    n = len(A)
    arcs = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    idx = {a: k for k, a in enumerate(arcs)}
    p = list(range(len(arcs)))
    def f(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    for g in auts:
        for (i, j) in arcs:
            a, b = f(idx[(i, j)]), f(idx[(g[i], g[j])])
            if a != b: p[a] = b
    return len({f(k) for k in range(len(arcs))})

def is_self_converse(A):
    n = len(A)
    outdeg = [sum(A[i]) for i in range(n)]
    indeg = [sum(A[j][i] for j in range(n)) for i in range(n)]
    out_by = {}
    for v in range(n):
        out_by.setdefault(outdeg[v], []).append(v)
    perm = [-1] * n; used = [False] * n
    order = sorted(range(n), key=lambda i: len(out_by.get(indeg[i], [])))
    def consistent(i, cand):
        for j in range(n):
            pj = perm[j]
            if pj == -1: continue
            if A[cand][pj] != A[j][i]: return False
            if A[pj][cand] != A[i][j]: return False
        return True
    def bt(k):
        if k == n: return True
        i = order[k]
        for cand in out_by.get(indeg[i], []):
            if used[cand]: continue
            if consistent(i, cand):
                perm[i] = cand; used[cand] = True
                if bt(k + 1): return True
                used[cand] = False; perm[i] = -1
        return False
    return bt(0)

def H_paths(A):
    n = len(A); full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            if not (mask >> v) & 1: continue
            cur = row[v]
            if cur == 0: continue
            Av = A[v]
            for w in range(n):
                if (mask >> w) & 1: continue
                if Av[w]: dp[mask | (1 << w)][w] += cur
    return sum(dp[full][v] for v in range(n))

def sym(name, A, H=None, compute_H=False):
    t0 = time.time()
    drt = is_DRT_adj(A)
    auts = automorphisms(A, cap=600000)
    vo = vertex_orbits(A, auts); ao = arc_orbits(A, auts)
    sc = is_self_converse(A)
    vt = (vo == 1); at = (vt and ao == 1)
    if compute_H and H is None:
        H = H_paths(A)
    det = det_I_plus_S(A)
    print(f"  {name:18s} DRT={'Y' if drt else 'n'} |Aut|={len(auts):8d} vO={vo:2d} arcO={ao:4d} "
          f"VT={'Y' if vt else 'n'} ArcT={'Y' if at else 'n'} SC={'Y' if sc else 'n'} "
          f"det(I+S)={det:>14}  H={H}  [{time.time()-t0:.1f}s]")
    return dict(name=name, drt=drt, naut=len(auts), vo=vo, ao=ao, sc=sc, vt=vt, at=at, H=H, det=det)

# ============================================================================
# n=7 : Paley=rotation? (at n=7, m=3, rotation S={1,2,3}, NOT Paley {1,2,4})
banner("n=7  rotation S={1,2,3} vs Paley S={1,2,4}")
A7r = circ(7, [1, 2, 3]); A7p = circ(7, [1, 2, 4])
sym("rotation7", A7r, compute_H=True)
sym("Paley7", A7p, compute_H=True)

# n=11
banner("n=11  rotation S={1..5} vs Paley QR_11")
A11r = circ(11, [1, 2, 3, 4, 5])
QR11 = sorted({(x * x) % 11 for x in range(1, 11)})
A11p = circ(11, QR11)
sym("rotation11", A11r, compute_H=True)
sym("Paley11", A11p, compute_H=True)

# n=15: rotation vs (no Paley). tower DRT done elsewhere; compute rotation H + det.
banner("n=15  rotation S={1..7}  (PALEY-FREE: no QR tournament; tower DRT is the DRT)")
A15r = circ(15, list(range(1, 8)))
print("  computing rotation15 symmetry + H (Held-Karp 2^15)...")
sym("rotation15", A15r, compute_H=True)

# n=19: rotation vs Paley QR_19. Use canon EXACT H (LEM-004); compute symmetry only.
banner("n=19  rotation S={1..9} vs Paley QR_19  (H from canon LEM-004; symmetry computed)")
A19r = circ(19, list(range(1, 10)))
QR19 = sorted({(x * x) % 19 for x in range(1, 19)})
A19p = circ(19, QR19)
H_rot19 = 1184212824763   # canon LEM-004 exact
H_pal19 = 1172695746915   # canon LEM-004 exact (Paley ranks 3rd)
sym("rotation19", A19r, H=H_rot19)
sym("Paley19", A19p, H=H_pal19)

banner("VERDICT SUMMARY")
print("""  CANON LEM-004: circulant H-max = Paley at n=7,11; = ROTATION at n=13,15,17,19.
  THM-472: det(I+S) Hadamard ceiling = DRT (spectral / weakly-regular) extremal.
  THM-499: H is STRICTLY FINER than the spectrum.

  Reading the symmetry columns above:
   * PALEY = arc-transitive (ArcT=Y), self-converse (SC=Y), a DRT, attains the
     det(I+S) ceiling. It is the MAXIMAL-SYMMETRY / SPECTRAL extremal.
   * ROTATION = the H-MAX (n>=13) but is NOT arc-transitive and NOT a DRT.
   * The half-arc-transitive (Doyle-Holt) object is the NS / 2-arc-orbit case:
     it is neither Paley nor rotation.

  => The arc-transitive (Doyle-Holt) lens characterizes the det(I+S)/DRT
     spectral ceiling, NOT the H-max. The two extremalities DECOUPLE at n>=13.
     Doyle-Holt is the OBSTRUCTION (wrong extremality), not the route, for H.""")

# ============================================================================
# APPENDIX (mac-mini-2026-06-21): the symmetry trichotomy of the DRT/H extremals
# across Paley-free orders n = 3 mod 4. Pins WHERE arc-transitivity is available.
# ============================================================================
def banner2(t): print("\n" + "=" * 74 + "\n  " + t + "\n" + "=" * 74)
from sympy import primefactors
banner2("APPENDIX: arc-transitivity availability across n = 3 mod 4")
print("  n   factors      prime-pow   nonab-order-n   symmetry available for the DRT")
def nonab(n):
    pf = primefactors(n)
    if len(pf) == 2 and n == pf[0]*pf[1]:
        p, q = sorted(pf); return (q-1) % p == 0
    return None
for n in [7,11,15,19,21,23,27,31,35,39,43,51,55]:
    if n % 4 != 3: continue
    pf = primefactors(n); pp = (len(pf) == 1)
    na = nonab(n)
    if pp:
        avail = "ARC-TRANSITIVE (Paley/QR circulant: |Aut|=n(n-1)/2, arcO=1)"
    elif na:
        avail = "HALF-ARC/NS only (non-ab carrier exists, e.g. F_21; VT but arcO>=2, NS)"
    else:
        avail = "DEFECTIVE only (no order-n group but Z_n; no circ DRT => non-VT NS)"
    print(f"  {n:3d} {str(pf):11s} {str(pp):>9s}   {str(na):>13s}   {avail}")
print("""
  TRICHOTOMY (proved/structural):
   * PRIME-POWER n (7,11,19,23,27,31,43): Paley/QR circulant DRT exists, is
     ARC-TRANSITIVE & self-converse. The det(I+S) Hadamard ceiling has a
     maximal-symmetry attainer. (But H-max is still rotation for n>=13.)
   * n = p*q with p | q-1 (21,39,55): no Paley, but a non-abelian Frobenius
     group of order n exists -> VERTEX-TRANSITIVE NON-SELF-CONVERSE tournament
     (kps HYP-2748, F_21). Half-arc-transitive regime: the converse Z_2 is
     UNREALIZED. Not arc-transitive.
   * n = p*q with p nmid q-1 (15,35,51): only Z_n; 0 circulant DRTs -> the DRT
     is FORCED non-VT, NS, DEFECTIVE symmetry (tower DRT(15): vO=3,arcO=7).
  So ARC-TRANSITIVITY (the Doyle-Holt 'symmetric' end) is available ONLY at
  prime-power orders. Everywhere else the DRT lives in the NS / half-arc-rigid
  or defective regime -- exactly the Doyle-Holt 'converse Z_2 unrealized' side.""")
