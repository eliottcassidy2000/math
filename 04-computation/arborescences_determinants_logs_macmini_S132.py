#!/usr/bin/env python3
"""
Arborescences on tournaments: the DETERMINANTAL RELAXATION of H, and its LOGARITHM
                                                        (mac-mini-2026-07-20-S132)
================================================================================
Owner: "merge in and extend ideas related to arborescences on tournaments and
determinants and logarithms."

MERGE NOTE.  The obvious merge -- "the characteristic polynomial and the OCF are the
same disjoint-cycle sum at different weights" -- is ALREADY REPO CANON: HYP-2514 /
THM-506's master cycle-packing polynomial Phi(T;{y_k}) = sum over linear subdigraphs
of prod y_{|C|}, with the spectrum its SIGNED vertex-graded (Sachs) face and H its
UNSIGNED ODD-ONLY face.  Not re-claimed.  What is genuinely absent from the repo is
ARBORESCENCES, so that is what this builds.

THE IDEA.  A Hamiltonian path from r IS a spanning out-arborescence rooted at r in
which every vertex has out-degree <= 1.  So arborescences are the RELAXATION of
Hamiltonian paths that drops the path constraint -- and unlike H, the relaxation is a
DETERMINANT (Tutte's matrix-tree theorem), hence polynomial time.  That is the exact
sense in which the determinant is the tractable shadow of the intractable H.

  # out-arborescences rooted at r  =  det of (D_in - A) with row/col r deleted
  sum_r a_r                        =  product of the NONZERO eigenvalues of D_in - A
                                       (Kirchhoff)  -- so the TOTAL is Laplacian-spectral

THE LOGARITHM.  log(sum_r a_r) = sum_i log(mu_i): the determinant becomes an additive
spectral sum.  Under ORDINAL SUM the two counts behave differently, and that difference
is the whole point:
  log H(T1 (+) T2)  = log H(T1) + log H(T2)                        -- no interaction
  log Sa(T1 (+) T2) = log Sa(T1) + log det(|T1|*I + L_in(T2))      -- SIZE-SHIFTED
"""
import numpy as np
from itertools import permutations, combinations
from fractions import Fraction

def scaffold(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return pairs, {p: k for k, p in enumerate(pairs)}, len(pairs)

def adj(code, pairs, n):
    A = np.zeros((n, n), dtype=np.int64)
    for e, (i, j) in enumerate(pairs):
        if code >> e & 1: A[j, i] = 1
        else:             A[i, j] = 1
    return A

def L_in(A):
    return np.diag(A.sum(axis=0)) - A          # D_in - A ;  1^T L_in = 0

def det_int(M):
    """exact integer determinant by fraction-free elimination."""
    M = [[Fraction(int(x)) for x in row] for row in M]
    n = len(M); det = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return 0
        if p != c: M[c], M[p] = M[p], M[c]; det = -det
        det *= M[c][c]
        inv = Fraction(1) / M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] * inv
            if f:
                for k in range(c, n): M[r][k] -= f * M[c][k]
    assert det.denominator == 1
    return int(det)

def arbor_matrixtree(A, r):
    L = L_in(A); idx = [i for i in range(len(A)) if i != r]
    return det_int(L[np.ix_(idx, idx)])

def arbor_bruteforce(A, r):
    """direct count: each v != r picks an in-neighbour parent; must form a tree rooted at r."""
    n = len(A); verts = [v for v in range(n) if v != r]
    parents = [[u for u in range(n) if A[u, v]] for v in range(n)]
    cnt = 0
    def rec(k, choice):
        nonlocal cnt
        if k == len(verts):
            # check acyclicity / reachability from r
            for v in verts:
                seen, x = set(), v
                while x != r:
                    if x in seen: return
                    seen.add(x); x = choice[x]
                    if x is None: return
            cnt += 1; return
        v = verts[k]
        for u in parents[v]:
            choice[v] = u; rec(k + 1, choice)
        choice[v] = None
    rec(0, {v: None for v in verts})
    return cnt

def ham_paths_from(A, r):
    n = len(A); dp = {(1 << r, r): 1}
    for _ in range(n - 1):
        nd = {}
        for (S, v), c in dp.items():
            for u in range(n):
                if S >> u & 1 or not A[v, u]: continue
                nd[(S | 1 << u, u)] = nd.get((S | 1 << u, u), 0) + c
        dp = nd
    return sum(dp.values())

def canon_codes(n):
    reps = {0}
    for k in range(2, n + 1):
        pk, ik, Ek = scaffold(k)
        op, _, _ = scaffold(k - 1)
        cand = []
        for r in reps:
            base = 0
            for e, (i, j) in enumerate(op):
                if r >> e & 1: base |= 1 << ik[(i, j)]
            for mask in range(1 << (k - 1)):
                v = base
                for b in range(k - 1):
                    if mask >> b & 1: v |= 1 << ik[(b, k - 1)]
                cand.append(v)
        p2 = (1 << np.arange(Ek, dtype=np.int64))
        Am = ((np.array(cand, dtype=np.int64)[:, None] >> np.arange(Ek)) & 1).astype(np.uint8)
        best = None
        for p in permutations(range(k)):
            src = np.empty(Ek, dtype=np.int64); fl = np.zeros(Ek, dtype=np.uint8)
            for e, (i, j) in enumerate(pk):
                a, b = p[i], p[j]
                t = ik[(min(a, b), max(a, b))]
                src[t] = e; fl[t] = 1 if a > b else 0
            c = (Am[:, src] ^ fl) @ p2
            best = c if best is None else np.minimum(best, c)
        reps = set(int(x) for x in best.tolist())
    return sorted(reps)

# ===================================================================== PART A
print("=" * 78)
print("PART A -- Tutte's matrix-tree on tournaments, and Kirchhoff's product")
print("=" * 78)
print(f"{'n':>3} {'classes':>8} {'matrixtree == brute force':>26} "
      f"{'sum_r a_r == prod nonzero mu':>29}")
for n in range(3, 7):
    pairs, idx, E = scaffold(n)
    reps = canon_codes(n)
    ok1 = ok2 = True
    for c in reps:
        A = adj(c, pairs, n)
        tot = 0
        for r in range(n):
            mt = arbor_matrixtree(A, r)
            if mt != arbor_bruteforce(A, r): ok1 = False
            tot += mt
        mu = np.linalg.eigvals(L_in(A).astype(float))
        nz = [m for m in mu if abs(m) > 1e-7]
        if abs(np.prod(nz).real - tot) > 1e-5 * max(1, tot): ok2 = False
    print(f"{n:>3} {len(reps):>8} {str(ok1):>26} {str(ok2):>29}")

print()
print("  Transitive tournament: L_in is upper triangular with in-degrees 0,1,...,n-1")
print("  on the diagonal, so sum_r a_r = (n-1)!  -- and ALL of it sits at the source,")
print("  since no arborescence can be rooted anywhere else.")
for n in range(3, 9):
    pairs, idx, E = scaffold(n)
    A = adj(0, pairs, n)
    tot = sum(arbor_matrixtree(A, r) for r in range(n))
    import math
    print(f"    n={n}: sum_r a_r = {tot:>8}   (n-1)! = {math.factorial(n-1):>8}   "
          f"H(TT_n) = {sum(ham_paths_from(A, r) for r in range(n))}")

print()
print("  Paley tournaments (q = 3 mod 4) are regular, so L_in = ((q-1)/2)I - A and")
print("  spec(A) = {(q-1)/2} u {(-1 +- i sqrt q)/2}, giving the CLOSED FORM")
print("      sum_r a_r = [ q(q+1)/4 ] ^ ((q-1)/2),     a_r = (1/q) * that.")
for q in (3, 7, 11, 19):
    QR = {(x * x) % q for x in range(1, q)}
    A = np.zeros((q, q), dtype=np.int64)
    for i in range(q):
        for k in QR: A[i, (i + k) % q] = 1
    tot = sum(arbor_matrixtree(A, r) for r in range(q))
    pred = (q * (q + 1) // 4) ** ((q - 1) // 2)
    H = sum(ham_paths_from(A, r) for r in range(q))
    print(f"    q={q:>2}: sum_r a_r = {tot:>18}  predicted {pred:>18}  match {tot==pred}"
          f"   H = {H}")

# ===================================================================== PART B
print()
print("=" * 78)
print("PART B -- arborescences RELAX Hamiltonian paths:  h_r <= a_r,  H <= sum_r a_r")
print("=" * 78)
print("  A Hamiltonian path from r IS a spanning out-arborescence rooted at r with every")
print("  out-degree <= 1.  Dropping that constraint is exactly the relaxation, and the")
print("  relaxation is a DETERMINANT while H is not (THM-505/506: only the SIGNED face of")
print("  the master cycle-packing polynomial collapses to det).")
print()
print(f"{'n':>3} {'h_r<=a_r always':>16} {'min ratio Sa/H':>16} {'max ratio Sa/H':>16} "
      f"{'argmin':>22} {'argmax':>22}")
for n in range(3, 8):
    pairs, idx, E = scaffold(n)
    reps = canon_codes(n)
    ok = True; best = None; worst = None
    for c in reps:
        A = adj(c, pairs, n)
        a = [arbor_matrixtree(A, r) for r in range(n)]
        h = [ham_paths_from(A, r) for r in range(n)]
        if any(hh > aa for hh, aa in zip(h, a)): ok = False
        ratio = sum(a) / sum(h)
        sc = sorted(A.sum(axis=1).tolist())
        if best is None or ratio < best[0]: best = (ratio, sc, sum(a), sum(h))
        if worst is None or ratio > worst[0]: worst = (ratio, sc, sum(a), sum(h))
    print(f"{n:>3} {str(ok):>16} {best[0]:>16.4f} {worst[0]:>16.4f} "
          f"{str(best[1]):>22} {str(worst[1]):>22}")
print("  (argmin/argmax reported as score sequences)")

# ===================================================================== PART C
print()
print("=" * 78)
print("PART C -- is sum_r a_r ADJACENCY-spectral?  (where it sits vs THM-499/500)")
print("=" * 78)
print("  Kirchhoff makes sum_r a_r LAPLACIAN-spectral by construction.  But L_in = D_in - A")
print("  knows the scores, and spec(A) pins only sum s_i and sum s_i^2 (via c3 = tr A^3/3).")
print("  So the question is whether ADJACENCY-cospectral tournaments can differ in sum a_r.")
print()
for n in range(4, 8):
    pairs, idx, E = scaffold(n)
    reps = canon_codes(n)
    byspec = {}
    for c in reps:
        A = adj(c, pairs, n)
        ch = tuple(np.round(np.poly(A.astype(float)), 6))
        byspec.setdefault(ch, []).append(c)
    split = 0; groups = 0
    for ch, cs in byspec.items():
        if len(cs) < 2: continue
        groups += 1
        vals = {sum(arbor_matrixtree(adj(c, pairs, n), r) for r in range(n)) for c in cs}
        if len(vals) > 1: split += 1
    print(f"  n={n}: {len(byspec)} adjacency-spectra over {len(reps)} classes; "
          f"{groups} cospectral groups; sum_r a_r DIFFERS inside {split} of them"
          f"  => adjacency-spectral? {split == 0 and groups > 0}")

# ===================================================================== PART D
print()
print("=" * 78)
print("PART D -- the LOGARITHM: additive with no interaction (H) vs size-shifted (arbor)")
print("=" * 78)
print("  Ordinal sum T1 (+) T2 makes L_in BLOCK UPPER TRIANGULAR:")
print("      [ L_in(T1)      -J   ]      so spec = spec(L1)  u  ( |T1| + spec(L2) )")
print("      [    0      pI + L2  ]")
print("  hence   sum_r a_r(T1 (+) T2) = sum_r a_r(T1) * det(|T1| I + L_in(T2)).")
print("  Taking logs:")
print("      log H (T1(+)T2) = log H(T1)  + log H(T2)                 -- NO interaction")
print("      log Sa(T1(+)T2) = log Sa(T1) + sum_mu log(|T1| + mu)     -- SIZE-SHIFTED")
print()
def ordinal(c1, n1, c2, n2):
    p1, i1, _ = scaffold(n1); p2, i2, _ = scaffold(n2)
    n = n1 + n2; pairs, idx, E = scaffold(n)
    A = np.zeros((n, n), dtype=np.int64)
    A1, A2 = adj(c1, p1, n1), adj(c2, p2, n2)
    A[:n1, :n1] = A1
    A[n1:, n1:] = A2
    A[:n1, n1:] = 1                                    # T1 beats all of T2
    return A

print(f"{'|T1|':>5} {'|T2|':>5} {'Sa(T)':>10} {'Sa(T1)*det(pI+L2)':>20} {'ok':>5} "
      f"{'H(T)':>7} {'H1*H2':>7} {'ok':>5}")
allok = True
for n1 in (2, 3):
    for n2 in (2, 3):
        p1, i1, e1 = scaffold(n1); p2, i2, e2 = scaffold(n2)
        for c1 in canon_codes(n1):
            for c2 in canon_codes(n2):
                A = ordinal(c1, n1, c2, n2)
                n = n1 + n2
                Sa = sum(arbor_matrixtree(A, r) for r in range(n))
                A1, A2 = adj(c1, p1, n1), adj(c2, p2, n2)
                Sa1 = sum(arbor_matrixtree(A1, r) for r in range(n1))
                pred = Sa1 * det_int(n1 * np.eye(n2, dtype=np.int64) + L_in(A2))
                H = sum(ham_paths_from(A, r) for r in range(n))
                H1 = sum(ham_paths_from(A1, r) for r in range(n1))
                H2 = sum(ham_paths_from(A2, r) for r in range(n2))
                ok1, ok2 = (Sa == pred), (H == H1 * H2)
                allok &= ok1 and ok2
                print(f"{n1:>5} {n2:>5} {Sa:>10} {pred:>20} {str(ok1):>5} "
                      f"{H:>7} {H1*H2:>7} {str(ok2):>5}")
print(f"  all ordinal-sum laws hold: {allok}")

print()
print("  TREE ENTROPY (1/n) log(sum_r a_r) -- the determinant made additive:")
print(f"{'n':>3} {'transitive':>12} {'max over classes':>18} {'min over classes':>18}")
import math
for n in range(3, 8):
    pairs, idx, E = scaffold(n)
    vals = []
    for c in canon_codes(n):
        A = adj(c, pairs, n)
        vals.append(sum(arbor_matrixtree(A, r) for r in range(n)))
    tt = sum(arbor_matrixtree(adj(0, pairs, n), r) for r in range(n))
    print(f"{n:>3} {math.log(tt)/n:>12.5f} {math.log(max(vals))/n:>18.5f} "
          f"{math.log(min(vals))/n:>18.5f}")

print()
print("SUMMARY")
print("  Arborescences are the DETERMINANTAL RELAXATION of Hamiltonian paths: same")
print("  spanning objects, out-degree constraint dropped, and the drop is exactly what")
print("  turns an intractable count into a determinant.  The logarithm then converts that")
print("  determinant into an additive spectral sum -- and the ordinal-sum laws show the")
print("  precise difference: H's log is additive with NO interaction term, the")
print("  arborescence log is additive with a SIZE-DEPENDENT SHIFT.")
