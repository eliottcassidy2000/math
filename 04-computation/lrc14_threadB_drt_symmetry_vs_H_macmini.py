#!/usr/bin/env python3
"""
lrc14_threadB_drt_symmetry_vs_H_macmini.py
mac-mini-2026-06-21 -- THREAD B (Doyle-Holt / Paley extremality), decisive test.

THE DECISIVE QUESTION. At n=7 (script _doyle_holt_paley_): among Z_7 circulants
ARC-TRANSITIVE <=> H-MAX, exactly (Paley pair, |Aut|=21=F_21). But n=7 has a
UNIQUE DRT (Paley). The real test needs n where:
  (i) MULTIPLE non-isomorphic DRTs exist (all cospectral, all attain the
      det(I+S) Hadamard ceiling THM-472), but
  (ii) H is STRICTLY FINER than the spectrum (THM-499), so different DRTs can
      have DIFFERENT H.
Then: does the ARC-TRANSITIVE DRT (Paley) have the LARGER H? Is half-arc /
NS (the Doyle-Holt structure) the LOWER-H sibling?

KNOWN (canon): the Mersenne-tower DRT at n=15 is NS, NOT vertex-transitive,
NOT arc-transitive (drt_mersenne_tower_kps1.out: VT=False, self-conv=False,
H=198335025). n=15 has NO Paley (15 not a prime power). So n=15 is the
PALEY-FREE DRT case: what is the H-max there, and is any DRT arc-transitive?

This script:
 (1) builds the Mersenne-tower DRT at n=7, 15 (THM-448) and QR-Paley DRT at
     n=7, 11, 19 (prime, 3 mod 4),
 (2) for EACH: |Aut|, vertex-orbits, arc-orbits, self-converse, arc-transitive,
     half-arc-analog (VT + 2 arc-orbits + NS), and H,
 (3) at n=11: confirms Paley QR_11 is the unique circulant DRT and reports its
     symmetry + H,
 (4) at n=15: builds the tower DRT AND searches all 2^7 circulant antisymmetric
     sets for any circulant DRT(15); compares H and symmetry. If the H-max DRT
     at n=15 is NS / not-arc-transitive, the Doyle-Holt 'arc-transitive => H-max'
     route FAILS in the Paley-free regime.
 (5) at n=19: QR_19 (Paley, arc-transitive) -- report its symmetry + H as the
     arc-transitive benchmark; note H computation cost.

VERDICT logic: Thread-B HYP (arc-transitive <=> H-max) holds iff at EVERY n the
unique H-max DRT is arc-transitive. A single NS/non-arc-transitive H-max DRT
(e.g. forced at n=15 where Paley is absent) is a clean OBSTRUCTION.
"""
import sys, time
from math import comb
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

def banner(t):
    print("\n" + "=" * 74 + "\n  " + t + "\n" + "=" * 74)

# ---------------- tournament utilities (adjacency 0/1) -----------------------
def is_tournament(A):
    n = len(A)
    for i in range(n):
        if A[i][i] != 0:
            return False
        for j in range(i + 1, n):
            if A[i][j] + A[j][i] != 1:
                return False
    return True

def skewM(A):
    n = len(A)
    return np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=np.int64)

def is_DRT_adj(A):
    """DRT iff S^2 = J - nI with S = A - A^T (equiv: every pair dominated by (n-3)/4)."""
    n = len(A)
    S = skewM(A).astype(np.int64)
    J = np.ones((n, n), dtype=np.int64)
    return np.array_equal(S @ S, J - n * np.eye(n, dtype=np.int64))

def det_I_plus_S(A):
    n = len(A)
    return round(float(np.linalg.det(np.eye(n) + skewM(A))))

def H_paths(A):
    n = len(A)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            cur = row[v]
            if cur == 0:
                continue
            Av = A[v]
            for w in range(n):
                if (mask >> w) & 1:
                    continue
                if Av[w]:
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[full][v] for v in range(n))

# ---------------- automorphisms (backtracking, degree pruning) ---------------
def automorphisms(A, cap=None):
    n = len(A)
    outdeg = [sum(A[i]) for i in range(n)]
    by_deg = {}
    for v in range(n):
        by_deg.setdefault(outdeg[v], []).append(v)
    perm = [-1] * n
    used = [False] * n
    order = sorted(range(n), key=lambda v: len(by_deg[outdeg[v]]))
    auts = []
    def consistent(i_orig, cand):
        for j_orig in range(n):
            pj = perm[j_orig]
            if pj == -1:
                continue
            if A[cand][pj] != A[i_orig][j_orig]:
                return False
            if A[pj][cand] != A[j_orig][i_orig]:
                return False
        return True
    def bt(k):
        if cap is not None and len(auts) >= cap:
            return
        if k == n:
            auts.append(perm[:]); return
        v = order[k]
        for cand in by_deg[outdeg[v]]:
            if used[cand]:
                continue
            if consistent(v, cand):
                perm[v] = cand; used[cand] = True
                bt(k + 1)
                used[cand] = False; perm[v] = -1
    bt(0)
    return auts

def vertex_orbits(A, auts):
    n = len(A); parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for p in auts:
        for v in range(n):
            ra, rb = find(v), find(p[v])
            if ra != rb: parent[ra] = rb
    return len({find(v) for v in range(n)})

def arc_orbits(A, auts):
    n = len(A)
    arcs = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    idx = {a: k for k, a in enumerate(arcs)}
    parent = list(range(len(arcs)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for p in auts:
        for (i, j) in arcs:
            ra, rb = find(idx[(i, j)]), find(idx[(p[i], p[j])])
            if ra != rb: parent[ra] = rb
    return len({find(k) for k in range(len(arcs))})

def is_self_converse(A):
    n = len(A)
    outdeg = [sum(A[i]) for i in range(n)]
    indeg = [sum(A[j][i] for j in range(n)) for i in range(n)]
    out_by = {}
    for v in range(n):
        out_by.setdefault(outdeg[v], []).append(v)
    perm = [-1] * n; used = [False] * n
    order = sorted(range(n), key=lambda i: len(out_by.get(indeg[i], [])))
    def consistent(i_orig, cand):
        for j_orig in range(n):
            pj = perm[j_orig]
            if pj == -1: continue
            if A[cand][pj] != A[j_orig][i_orig]: return False
            if A[pj][cand] != A[i_orig][j_orig]: return False
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

def report(name, A, do_H=True, cap_aut=400000):
    t0 = time.time()
    drt = is_DRT_adj(A)
    auts = automorphisms(A, cap=cap_aut)
    vo = vertex_orbits(A, auts); ao = arc_orbits(A, auts)
    sc = is_self_converse(A)
    vt = (vo == 1); at = (vt and ao == 1); half = (vt and ao == 2 and not sc)
    H = H_paths(A) if do_H else None
    det = det_I_plus_S(A)
    print(f"  {name:26s} DRT={'Y' if drt else 'n'}  |Aut|={len(auts):7d} vO={vo:2d} arcO={ao:3d} "
          f"VT={'Y' if vt else 'n'} ArcT={'Y' if at else 'n'} SelfConv={'Y' if sc else 'n'} "
          f"halfArc={'Y' if half else 'n'}  H={H}  det(I+S)={det}  [{time.time()-t0:.1f}s]")
    return dict(name=name, drt=drt, naut=len(auts), vo=vo, ao=ao, sc=sc, vt=vt,
                at=at, half=half, H=H, det=det)

# ---------------- constructions ---------------------------------------------
def paley_qr(p):
    """Paley/QR circulant tournament on Z_p (p prime, 3 mod 4). A[i][j]=1 iff j-i in QR."""
    QR = set((x * x) % p for x in range(1, p))
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in QR:
                A[i][j] = 1
    return A

def mersenne_tower_drt(k):
    """Core DRT on n=2^k-1 from the skew-doubling tower S -> [[S,S],[S-2I,2I-S]],
    seed S2 = [[1,1],[-1,1]]. Core tournament = strip border, A = (M-... )."""
    S = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(k - 1):
        n = S.shape[0]; I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
    # S is skew-Hadamard of order 2^k. Normalize first row to all +1, then the
    # core on rows/cols 1..(2^k-1) gives the DRT skew matrix S0; A from sign.
    # Normalize: multiply col j and row j by S[0,j] so first row/col are +1.
    m = S.shape[0]
    d = np.array([S[0, j] for j in range(m)], dtype=np.int64)
    Sn = (d[:, None] * S * d[None, :])
    # now Sn[0,:] = +1, Sn[:,0]=+1 (since seed gives that after normalization)
    core = Sn[1:, 1:]  # (2^k-1)x(2^k-1) skew matrix with +-1 off diagonal
    n = core.shape[0]
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and core[i, j] == 1:
                A[i][j] = 1
            elif i != j and core[i, j] == -1:
                A[i][j] = 0
    # ensure tournament
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] + A[j][i] != 1:
                # fix from skew sign
                A[i][j] = 1 if core[i, j] > 0 else 0
                A[j][i] = 1 - A[i][j]
    return A

def all_circulant_drts(n):
    """All circulant tournaments on Z_n that are DRTs (n odd). Antisymmetric S:
    one of each {s, n-s}. n must be such that S is a skew-Hadamard difference set."""
    m = (n - 1) // 2
    found = []
    halves = list(range(1, m + 1))  # representatives 1..m; partner n-s
    for bits in range(1 << m):
        S = []
        for k in range(m):
            s = halves[k]
            S.append(s if (bits >> k) & 1 else (n - s))
        A = [[0] * n for _ in range(n)]
        Sset = set(S)
        for i in range(n):
            for j in range(n):
                if i != j and (j - i) % n in Sset:
                    A[i][j] = 1
        if is_tournament(A) and is_DRT_adj(A):
            found.append((tuple(sorted(S)), A))
    return found

# ============================================================================
banner("(B) n=7 -- unique DRT (Paley QR_7) vs Mersenne-tower DRT(7): same object?")
A7p = paley_qr(7)
A7t = mersenne_tower_drt(3)  # 2^3-1 = 7
r7p = report("Paley QR_7", A7p)
r7t = report("Tower DRT(7)", A7t)

banner("(C) n=11 -- Paley QR_11 (unique circulant DRT)")
A11 = paley_qr(11)
r11 = report("Paley QR_11", A11)
cd11 = all_circulant_drts(11)
print(f"  #circulant DRTs at n=11 (over antisymmetric sets): {len(cd11)}")
for Sset, _ in cd11:
    print(f"    circulant DRT set S={list(Sset)}")

banner("(D) n=15 -- THE PALEY-FREE CASE: tower DRT (NS) + any circulant DRT")
A15t = mersenne_tower_drt(4)  # 2^4-1 = 15
print("  building tower DRT(15) symmetry + H (may take a few seconds)...")
r15t = report("Tower DRT(15)", A15t)
print("  searching all 2^7=128 antisymmetric circulant sets on Z_15 for DRTs...")
cd15 = all_circulant_drts(15)
print(f"  #circulant DRTs at n=15: {len(cd15)}")
r15_circ = []
for Sset, Ac in cd15:
    rr = report(f"circ DRT15 S={list(Sset)}", Ac)
    r15_circ.append(rr)

# Collect all known DRTs at n=15 and find the H-max
all15 = [r15t] + r15_circ
if all15 and all(r['H'] is not None for r in all15):
    Hmax15 = max(r['H'] for r in all15)
    print(f"\n  H-max among KNOWN DRTs at n=15 = {Hmax15}")
    for r in all15:
        if r['H'] == Hmax15:
            print(f"    H-MAX DRT: {r['name']}  ArcTrans={r['at']} SelfConv={r['sc']} VT={r['vt']}")
    print("  >> Is the H-max DRT at n=15 ARC-TRANSITIVE?",
          any(r['at'] for r in all15 if r['H'] == Hmax15))
    print("  >> Is the H-max DRT at n=15 SELF-CONVERSE (SC)?",
          any(r['sc'] for r in all15 if r['H'] == Hmax15))

banner("SYNTHESIS TABLE")
print(f"  {'object':26s} {'n':>3} {'DRT':>4} {'VT':>3} {'ArcT':>5} {'SC':>3} {'H':>14}")
for r in [r7p, r7t, r11, r15t] + r15_circ:
    if r is None: continue
    nn = {'Paley QR_7':7,'Tower DRT(7)':7,'Paley QR_11':11,'Tower DRT(15)':15}.get(r['name'],15)
    print(f"  {r['name']:26s} {nn:>3} {('Y' if r['drt'] else 'n'):>4} "
          f"{('Y' if r['vt'] else 'n'):>3} {('Y' if r['at'] else 'n'):>5} "
          f"{('Y' if r['sc'] else 'n'):>3} {str(r['H']):>14}")
