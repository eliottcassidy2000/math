#!/usr/bin/env python3
"""
Vapnik-Chervonenkis dimension on both sides of the odd/even axis  (mac-mini-S131)
================================================================================
Owner: "think vapnik-chervonenkis and chase the high leverage question, see the
relation between odd valued functions and tournament adjacent ideas.  they both
relate also to even concepts like even graphs and even functions."

The high-leverage question turns out to be: WHY is VC trivial on one side of the
repo's odd/even axis and hard on the other?  Answer tested here: LINEARITY.

  EVEN SIDE (symmetric functions, F_2-LINEAR).  Even graphs = the cycle space; the
  switching group = the cut space.  For ANY F_2-linear code C, a coordinate set S is
  shattered iff the generator columns on S are independent, so
        VC(C) = dim(C)   and   shattered sets = INDEPENDENT SETS OF A MATROID.
  For the cut space of K_n the matroid is graphic: the shattered sets are the FORESTS.
  VC carries no information beyond the dimension.  And a coset behaves identically,
  so EVERY SWITCHING CLASS has VC = n-1 with the forests as its shattered sets.

  ODD SIDE (skew functions, NOT linear).  The out-neighbourhood system {N+(v)} of a
  single tournament.  Shattering a k-set needs 2^k vertices, so VC <= floor(log2 n);
  the bound is attained.  Here VC is a genuine combinatorial quantity, in Erdos /
  Schutte territory -- and by THM-1420 the tournament side admits NO F_2-linear
  invariants at all, so it CANNOT be linearised into the trivial case.
"""
import numpy as np
from itertools import permutations, combinations, chain

def scaffold(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return pairs, {p: k for k, p in enumerate(pairs)}, len(pairs)

def subsets(it):
    s = list(it)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

# ------------------------------------------------------------------ linear side
def rank_f2(cols):
    """rank over F_2 of a list of integer bitmask columns."""
    basis = []
    for v in cols:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v); basis.sort(reverse=True)
    return len(basis)

def shattered_by_code(gen, N):
    """gen = list of generator rows (bitmasks over N coords).  S is shattered iff the
       projection of the code onto S is all of F_2^S iff the generator COLUMNS on S
       are F_2-independent."""
    cols = []
    for c in range(N):
        v = 0
        for r, g in enumerate(gen):
            if g >> c & 1: v |= 1 << r
        cols.append(v)
    out = []
    for S in subsets(range(N)):
        if rank_f2([cols[c] for c in S]) == len(S): out.append(S)
    return out, cols

def brute_shatter_code(gen, N):
    """independent check: enumerate the code and test shattering directly."""
    code = set()
    for m in range(1 << len(gen)):
        v = 0
        for b in range(len(gen)):
            if m >> b & 1: v ^= gen[b]
        code.add(v)
    out = []
    for S in subsets(range(N)):
        pats = {tuple((w >> c) & 1 for c in S) for w in code}
        if len(pats) == 1 << len(S): out.append(S)
    return out

print("=" * 78)
print("PART A -- on the EVEN/LINEAR side, VC is just the dimension, and the shattered")
print("          sets are the INDEPENDENT SETS OF A MATROID.  For the cut space of K_n")
print("          (= the switching group) that matroid is graphic: shattered = FORESTS.")
print("=" * 78)
print(f"{'n':>3} {'edges':>6} {'cut dim':>8} {'VC(cut)':>8} {'shattered = forests?':>21} "
      f"{'cycle dim':>10} {'VC(cycle)':>10} {'brute agrees':>13}")
for n in range(3, 7):
    pairs, idx, E = scaffold(n)
    stars = []
    for v in range(n):
        b = 0
        for e, (i, j) in enumerate(pairs):
            if i == v or j == v: b |= 1 << e
        stars.append(b)
    sh_cut, _ = shattered_by_code(stars, E)
    # cycle space generators: fundamental cycles w.r.t. a spanning star at vertex 0
    tree = {idx[(0, v)] for v in range(1, n)}
    cyc = []
    for e, (i, j) in enumerate(pairs):
        if e in tree: continue
        cyc.append((1 << e) | (1 << idx[(0, i)]) | (1 << idx[(0, j)]))
    sh_cyc, _ = shattered_by_code(cyc, E)
    # is "shattered" exactly "forest"?
    def is_forest(S):
        par = list(range(n))
        def f(a):
            while par[a] != a: par[a] = par[par[a]]; a = par[a]
            return a
        for c in S:
            i, j = pairs[c]
            ri, rj = f(i), f(j)
            if ri == rj: return False
            par[ri] = rj
        return True
    forests = [S for S in subsets(range(E)) if is_forest(S)]
    agree = set(sh_cut) == set(forests)
    brute_ok = (set(brute_shatter_code(stars, E)) == set(sh_cut) and
                set(brute_shatter_code(cyc, E)) == set(sh_cyc))
    vcc = max(len(S) for S in sh_cut); vcz = max(len(S) for S in sh_cyc)
    print(f"{n:>3} {E:>6} {n-1:>8} {vcc:>8} {str(agree):>21} "
          f"{E-n+1:>10} {vcz:>10} {str(brute_ok):>13}")

print()
print("  So every SWITCHING CLASS (a coset of the cut space) has VC = n-1 exactly,")
print("  with the forests of K_n as its shattered sets -- translation cannot change")
print("  surjectivity of a projection, so cosets and the code agree.")

# ------------------------------------------------------------------ odd side
def out_nbrs(code, pairs, n):
    N = [0] * n
    for e, (i, j) in enumerate(pairs):
        if code >> e & 1: N[j] |= 1 << i
        else:             N[i] |= 1 << j
    return N

def vc_tournament(code, pairs, n):
    N = out_nbrs(code, pairs, n)
    best = 0
    for k in range(1, n + 1):
        if (1 << k) > n: break
        found = False
        for S in combinations(range(n), k):
            pats = {tuple((N[v] >> s) & 1 for s in S) for v in range(n)}
            if len(pats) == (1 << k): found = True; break
        if found: best = k
        else: break
    return best

def schutte_k(code, pairs, n):
    """largest k such that EVERY k-set is beaten by some vertex (property S_k)."""
    N = out_nbrs(code, pairs, n)
    best = 0
    for k in range(1, n):
        ok = True
        for S in combinations(range(n), k):
            m = 0
            for s in S: m |= 1 << s
            if not any((N[v] & m) == m for v in range(n)): ok = False; break
        if ok: best = k
        else: break
    return best

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
        A = ((np.array(cand, dtype=np.int64)[:, None] >> np.arange(Ek)) & 1).astype(np.uint8)
        best = None
        for p in permutations(range(k)):
            src = np.empty(Ek, dtype=np.int64); fl = np.zeros(Ek, dtype=np.uint8)
            for e, (i, j) in enumerate(pk):
                a, b = p[i], p[j]
                t = ik[(min(a, b), max(a, b))]
                src[t] = e; fl[t] = 1 if a > b else 0
            c = (A[:, src] ^ fl) @ p2
            best = c if best is None else np.minimum(best, c)
        reps = set(int(x) for x in best.tolist())
    return sorted(reps)

print()
print("=" * 78)
print("PART B -- on the ODD/NONLINEAR side, VC is a real combinatorial quantity")
print("=" * 78)
print(f"{'n':>3} {'classes':>8} {'floor(log2 n)':>14} {'max VC':>7} {'VC distribution':>28} "
      f"{'max Schutte k':>14}")
for n in range(3, 8):
    pairs, idx, E = scaffold(n)
    reps = canon_codes(n)
    dist = {}
    sch = 0
    for r in reps:
        v = vc_tournament(r, pairs, n); dist[v] = dist.get(v, 0) + 1
        sch = max(sch, schutte_k(r, pairs, n))
    import math
    lg = int(math.log2(n))
    print(f"{n:>3} {len(reps):>8} {lg:>14} {max(dist):>7} "
          f"{str(sorted(dist.items())):>28} {sch:>14}")

print()
print("  Upper bound is immediate: shattering a k-set needs 2^k distinct patterns and")
print("  there are only n vertices, so VC <= floor(log2 n).")
print("  Lower bound by construction: put a TRANSITIVE tournament on S = {0..k-1}, whose")
print("  in-S out-neighbourhoods are the k distinct nested sets {1..k-1} ... {} ; assign")
print("  the remaining 2^k - k patterns freely to outside vertices.  So max VC = floor(log2 n).")

# explicit n=8 witness for VC = 3
n = 8; pairs, idx, E = scaffold(n)
S = (0, 1, 2)
want = [tuple((m >> b) & 1 for b in range(3)) for m in range(8)]
code = 0
for e, (i, j) in enumerate(pairs):
    if i in S and j in S: continue                 # S transitive: i->j for i<j, bit 0
for e, (i, j) in enumerate(pairs):
    if i in S and j not in S:
        # vertex j (outside) gets pattern index j (3..7) plus we need 3 more
        pass
# build directly: vertex v outside gets the pattern of index v
outer = [v for v in range(n) if v not in S]
insid = {0: (0, 1, 1), 1: (0, 0, 1), 2: (0, 0, 0)}   # transitive 0->1->2, 0->2
assign = dict(insid)
used = set(assign.values())
free = [w for w in want if w not in used]
for v, w in zip(outer, free): assign[v] = w
for e, (i, j) in enumerate(pairs):
    if i in S and j in S:
        if not (i < j): code |= 1 << e             # transitive, already i->j
    elif i in S and j not in S:
        if assign[j][i] == 1: code |= 1 << e       # j beats i  => bit set
    elif j in S and i not in S:
        if assign[i][j] == 0: code |= 1 << e       # j beats i
print()
print(f"  explicit n=8 witness: VC = {vc_tournament(code, pairs, n)} "
      f"(floor(log2 8) = 3), code = {code}")

print()
print("=" * 78)
print("PART C -- where does VC sit in the invariant tiers?")
print("=" * 78)
for n in range(4, 8):
    pairs, idx, E = scaffold(n)
    reps = canon_codes(n)
    switches = []
    for m in range(1 << (n - 1)):
        b = 0
        for e, (i, j) in enumerate(pairs):
            bi = (m >> i & 1) if i < n - 1 else 0
            bj = (m >> j & 1) if j < n - 1 else 0
            if bi ^ bj: b |= 1 << e
        switches.append(b)
    moved = 0
    for r in reps:
        v0 = vc_tournament(r, pairs, n)
        if any(vc_tournament(r ^ s, pairs, n) != v0 for s in switches): moved += 1
    print(f"  n={n}: VC changes under switching for {moved}/{len(reps)} classes "
          f"=> switching-invariant? {moved == 0}")

print()
print("SUMMARY -- the high-leverage answer")
print("  LINEARITY TRIVIALISES VC.  On the even/linear side (cut space = switching group,")
print("  cycle space = even graphs) VC equals the dimension and the shattered sets are")
print("  matroid independent sets -- forests, for the switching group.  Nothing new is")
print("  measured.  On the odd/nonlinear side VC is bounded by floor(log2 n), attains it,")
print("  and varies across tournaments.  THM-1420 says the odd side admits NO F_2-linear")
print("  invariants, so it can never be pushed into the trivial linear case.")
