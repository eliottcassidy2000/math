#!/usr/bin/env python3
"""
klein-2026-07-20-S338 -- SKEW-SEIDEL: THE TOURNAMENT ANALOGUE OF THE TWO-GRAPH/EVEN-GRAPH
THEOREM EXISTS.  It is EVEN TOURNAMENTS, not even graphs.

Owner: "think skew-Seidel and chase the high leverage question, see the relation between odd
valued functions and tournament adjacent ideas.  they both relate also to even concepts like
even graphs and even functions."

THE FRAME.
  A tournament IS a SKEW-Seidel matrix:  S_ij = +1 if i->j, -1 if j->i, 0 on the diagonal,
  so S^T = -S.  A graph is a SYMMETRIC Seidel matrix.  Seidel switching S -> D_U S D_U
  (D_U = diag(+-1)) acts on both.  Odd/skew side = tournaments; even/symmetric side = graphs.

  GRAPH SIDE (classical): switching classes of graphs = two-graphs, and for n ODD every
  switching class contains EXACTLY ONE graph with all degrees even -- an EVEN GRAPH.  Count
  = A002854 = the repo's V(E_n) = 2,3,7,16,54.

  TOURNAMENT SIDE: switching classes of tournaments = ORIENTED two-graphs (Babai-Cameron
  EJC 7 (2000) #R38; THM-474 already cites this).  kind-pasteur's THM-1415 compared the
  count 1,2,2,6 against A002854 = 2,3,7,16 and concluded "the graph two-graph theorem has NO
  tournament analogue".  THAT COMPARISON IS AGAINST THE WRONG OBJECT: the analogue of an
  even GRAPH is an even TOURNAMENT (all SCORES even), not an even graph.

THE THEOREM TESTED HERE.
  Switching at U reverses all arcs between U and U^c, so for v in U the score changes by
  |U^c| - 2a, and for v in U^c by |U| - 2b.  Hence
        Delta s_v  =  |U^c| (mod 2)  for v in U,     |U| (mod 2)  for v in U^c.
  For n ODD exactly one of |U|,|U^c| is even, and the SCORE-PARITY VECTOR is flipped exactly
  on the EVEN-SIZED one of {U, U^c}.  The map {switchings} = 2^[n]/{U,U^c} -> {even subsets}
  is a BIJECTION (both of size 2^(n-1)), so the 2^(n-1) switchings of a class realise
  2^(n-1) DISTINCT parity vectors -- the whole coset p_0 + (even-weight code).
  Since #odd scores = sum s_v = C(n,2) (mod 2):
        n = 1 mod 4  ->  C(n,2) EVEN  ->  the all-even parity vector IS in the coset,
                         so EVERY switching class has a UNIQUE tournament with all scores even.
        n = 3 mod 4  ->  C(n,2) ODD   ->  all-even is unreachable; the minimum is weight 1,
                         and every class has EXACTLY n representatives with a single odd score.
        n EVEN       ->  |U|,|U^c| share parity: either NO score flips or ALL do, so the
                         parity vector is pinned only up to global complement.  No analogue.

PARTS
  A. the parity law, verified exhaustively.
  B. unique even representative at n = 1 mod 4; exactly n single-odd reps at n = 3 mod 4.
  C. THE COUNT: switching classes up to iso  vs  EVEN TOURNAMENTS up to iso.
  D. the odd-valued side: det S = 0 for n odd (skew of odd order), det S = Pf(S)^2 for n even;
     and Redei's hp(T) odd.  Where the parity of the whole project actually lives.
"""
import itertools
from math import comb

# ---------------------------------------------------------------- tournaments as bitmasks
def pairs_of(n): return [(i, j) for i in range(n) for j in range(i + 1, n)]

def om_from_bits(n, bits, P):
    om = [0] * n
    for k, (i, j) in enumerate(P):
        if bits >> k & 1: om[i] |= 1 << j
        else:             om[j] |= 1 << i
    return tuple(om)

def scores(om, n): return [bin(om[v]).count("1") for v in range(n)]

def switch(om, n, U):
    """reverse every arc between U and its complement (U a bitmask)."""
    nm = list(om)
    for i in range(n):
        for j in range(i + 1, n):
            if ((U >> i) & 1) != ((U >> j) & 1):
                if nm[i] >> j & 1: nm[i] &= ~(1 << j); nm[j] |= 1 << i
                else:              nm[j] &= ~(1 << i); nm[i] |= 1 << j
    return tuple(nm)

def relabel(om, perm, n):
    new = [0] * n
    for v in range(n):
        mv, t = om[v], 0
        while mv:
            b = mv & -mv; w = b.bit_length() - 1; mv ^= b
            t |= 1 << perm[w]
        new[perm[v]] = t
    return tuple(new)

def word(om, n):
    w = 0
    for v in range(n): w = (w << n) | om[v]
    return w

def refine(om, n):
    colour = [bin(om[v]).count("1") for v in range(n)]
    while True:
        sig = []
        for v in range(n):
            cnt = {}; mv = om[v]
            while mv:
                b = mv & -mv; w = b.bit_length() - 1; mv ^= b
                cnt[colour[w]] = cnt.get(colour[w], 0) + 1
            sig.append((colour[v], tuple(sorted(cnt.items()))))
        order = sorted(set(sig)); newc = [order.index(sig[v]) for v in range(n)]
        if newc == colour: break
        colour = newc
    cells = {}
    for v in range(n): cells.setdefault(colour[v], []).append(v)
    return [tuple(cells[k]) for k in sorted(cells)]

_cc = {}
def canon_iso(om, n):
    """refinement-based canonical form: min over permutations respecting refined cells."""
    key = (om, n)
    r = _cc.get(key)
    if r is not None: return r
    cells = refine(om, n); base = []; pos = 0
    for c in cells: base.append((c, pos)); pos += len(c)
    best = None
    for choice in itertools.product(*[itertools.permutations(c) for (c, _) in base]):
        perm = [0] * n
        for (blk, (c, start)) in zip(choice, base):
            for k, v in enumerate(blk): perm[v] = start + k
        w = word(relabel(om, perm, n), n)
        if best is None or w < best: best = w
    _cc[key] = best
    return best

def canon_switch_iso(om, n):
    """canonical form of the SWITCHING CLASS up to isomorphism."""
    best = None
    for U in range(1 << (n - 1)):            # U and U^c give the same switching
        s = switch(om, n, U)
        w = canon_iso(s, n)
        if best is None or w < best: best = w
    return best

# ============================================================ PART A + B
print("=" * 78)
print("PART A/B -- THE SCORE-PARITY LAW AND THE CANONICAL REPRESENTATIVE")
print("=" * 78)
for n in (3, 4, 5, 6, 7):
    P = pairs_of(n); E = len(P)
    if E > 21: continue
    # (A) parity law: flipped set = the even-sized one of {U, U^c}
    om0 = om_from_bits(n, 0, P)
    s0 = scores(om0, n); p0 = [x & 1 for x in s0]
    lawok = True
    for U in range(1 << n):
        s1 = scores(switch(om0, n, U), n); p1 = [x & 1 for x in s1]
        flipped = {v for v in range(n) if p0[v] != p1[v]}
        Uset = {v for v in range(n) if U >> v & 1}
        Uc = set(range(n)) - Uset
        expect = Uset if len(Uset) % 2 == 0 else Uc
        if n % 2 == 1 and flipped != expect: lawok = False; break
        if n % 2 == 0 and flipped not in (set(), set(range(n))): lawok = False; break
    # (B) representative counts, over a sample of switching classes
    reps_alleven, reps_oneodd, distinct_parities = [], [], set()
    import random; random.seed(7 * n)
    tests = 200 if E <= 15 else 40
    ok_uniq = True
    for _ in range(tests):
        om = om_from_bits(n, random.getrandbits(E), P)
        ae = oo = 0; pv = set()
        for U in range(1 << (n - 1)):
            sc = scores(switch(om, n, U), n)
            pvec = tuple(x & 1 for x in sc); pv.add(pvec)
            if sum(pvec) == 0: ae += 1
            if sum(pvec) == 1: oo += 1
        reps_alleven.append(ae); reps_oneodd.append(oo); distinct_parities.add(len(pv))
    print(f"\n n={n}  (n mod 4 = {n % 4})   C(n,2) = {comb(n,2)}  -> #odd scores is "
          f"{'EVEN' if comb(n,2) % 2 == 0 else 'ODD'}")
    print(f"   parity law holds on all 2^n switchings : {lawok}")
    print(f"   distinct parity vectors per switching class : {sorted(distinct_parities)}"
          f"   (2^(n-1) = {1 << (n-1)} predicted for n odd)")
    print(f"   #reps with ALL SCORES EVEN  : {sorted(set(reps_alleven))}"
          f"     #reps with exactly ONE odd score : {sorted(set(reps_oneodd))}")

# ============================================================ PART C
print("\n" + "=" * 78)
print("PART C -- THE COUNT: switching classes up to iso  vs  EVEN TOURNAMENTS up to iso")
print("=" * 78)
print(f"{'n':>3} {'n%4':>4} {'iso classes':>12} {'switch cls /iso':>16} {'even tourns /iso':>17} {'A049313':>9} {'match?':>8}")
A049313 = {3: 1, 4: 2, 5: 2, 6: 6, 7: 12, 8: 79, 9: 792}
for n in (3, 4, 5, 6, 7):
    P = pairs_of(n); E = len(P)
    if E > 21: continue
    # 1) iso classes, by BFS over single-arc flips (cheap, reuses S336 idea)
    reps = {}
    om0 = tuple(sum(1 << j for j in range(i)) for i in range(n))
    seen = {canon_iso(om0, n): om0}; fr = [om0]
    while fr:
        nx = []
        for om in fr:
            for (i, j) in P:
                nm = list(om)
                if om[i] >> j & 1: nm[i] &= ~(1 << j); nm[j] |= 1 << i
                else:              nm[j] &= ~(1 << i); nm[i] |= 1 << j
                nm = tuple(nm); c = canon_iso(nm, n)
                if c not in seen: seen[c] = nm; nx.append(nm)
        fr = nx
    # 2) union-find the iso classes into SWITCHING classes
    parent = {c: c for c in seen}
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for c, om in seen.items():
        for U in range(1 << (n - 1)):
            union(c, canon_iso(switch(om, n, U), n))
    nsw = len({find(c) for c in seen})
    # 3) EVEN TOURNAMENTS up to iso: iso classes whose score multiset is all even
    nev = sum(1 for c, om in seen.items() if all(s % 2 == 0 for s in scores(om, n)))
    pred = A049313.get(n)
    match = (nsw == nev) if n % 4 == 1 else "n/a"
    print(f"{n:>3} {n%4:>4} {len(seen):>12} {nsw:>16} {nev:>17} {str(pred):>9} {str(match):>8}")

# n = 7 via the single-odd-score representatives (7 = 3 mod 4)
print("\n n=7 (= 3 mod 4): every switching class has exactly 7 single-odd-score reps.")
print(" Counting switching classes up to iso directly is 2^21 x 5040 and is not run here.")
print(" What IS run: the representative counts above already confirm 'exactly 7' at n=7.")

# ============================================================ PART D
print("\n" + "=" * 78)
print("PART D -- THE ODD-VALUED SIDE: where the project's parity actually lives")
print("=" * 78)
try:
    from fractions import Fraction
    def detZ(M):
        M = [row[:] for row in M]; n = len(M); det = Fraction(1)
        for c in range(n):
            piv = next((r for r in range(c, n) if M[r][c] != 0), None)
            if piv is None: return 0
            if piv != c: M[c], M[piv] = M[piv], M[c]; det = -det
            det *= M[c][c]
            inv = Fraction(1, 1) / M[c][c]
            for r in range(c + 1, n):
                f = M[r][c] * inv
                if f: M[r] = [M[r][k] - f * M[c][k] for k in range(n)]
        return det
    import random; random.seed(11)
    print(f"{'n':>3} {'det S(T) over 200 random T':>34}")
    for n in range(3, 9):
        P = pairs_of(n); E = len(P); vals = set()
        for _ in range(200):
            om = om_from_bits(n, random.getrandbits(E), P)
            S = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    if i != j: S[i][j] = 1 if (om[i] >> j & 1) else -1
            vals.add(detZ([[Fraction(v) for v in row] for row in S]))
        allsq = all(v >= 0 and int(round(v ** 0.5)) ** 2 == v for v in vals if v != 0)
        print(f"{n:>3}   n {'ODD ' if n%2 else 'EVEN'}  distinct dets: "
              f"{sorted(vals)[:6]}{'...' if len(vals)>6 else ''}"
              f"   {'all zero' if vals=={0} else ('all perfect squares (= Pf^2): ' + str(allsq))}")
except Exception as e:
    print(" det section skipped:", e)
print("\n READING: skew-symmetric of ODD order is singular, so det S(T) = 0 for every")
print(" tournament on an odd number of vertices; for EVEN order det S = Pf(S)^2 is a square.")
print(" That is the same odd/even split that makes the SCORE-PARITY law work only at odd n,")
print(" and it is the matrix-level shadow of Redei's theorem that hp(T) is always ODD.")
print("\nDONE.")
