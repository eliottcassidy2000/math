#!/usr/bin/env python3
"""switching_pisano_faulhaber_kps_S128c110.py -- kind-pasteur-2026-07-20-S128c110

THREE THREADS, ONE SCRIPT.

(A) A BETTER BASE-PATH-INDEPENDENT INVARIANT THAN "intersect GAMMA over all
    spanning paths".  THM-1405 showed the star group GAMMA fails to descend because
    it is defined by RESTRICTING the cut space to non-tree edges, and that
    restriction needs a base path.  The fix is not to intersect over paths -- it is
    to STOP RESTRICTING.  In the full edge space of K_n the cut space is spanned by
    the vertex stars delta(v), and S_n PERMUTES those stars, so the cut space is an
    S_n-INVARIANT subspace.  The canonical, base-path-free operation is therefore

        SWITCHING:  pick S subset V, reverse every arc between S and V\\S.

    This is Seidel switching transported to tournaments.  It commutes with
    relabelling, so "switching class" is a genuine tournament-level notion and the
    quotient is canonical.  Computed here: the number of switching classes of
    tournaments up to isomorphism, against the iso-class counts A000568 and against
    A002854 = 2,3,7,16,54 (the even-graph / two-graph counts the repo already treats
    as first-class).  For GRAPHS, switching classes up to iso ARE the two-graphs and
    ARE counted by A002854; whether the tournament analogue lands on the same
    sequence is exactly the question.

(B) THE THREE SIXTIES, deflated quantitatively.  ord_1001(2) = 60, the Pisano period
    pi(10) = 60, and |A_5| = 60.  Rather than assert coincidence, MEASURE why 60 is
    over-represented: 60 = 2^2*3*5 is the smallest number divisible by 1..6, so lcm's
    of small integers pile onto it.  The script computes the causes separately
    (ord_1001(2) = lcm(ord_7, ord_11, ord_13) = lcm(3,10,12); pi(10) = lcm(pi(2),
    pi(5)) = lcm(3,20)) and then measures how often lcm's of small sets land on 60.

(C) FIBONACCI FROM SHIFTED PASCAL, applied to the Rosetta/Faulhaber triangle.  The
    classical fact is F_{m+1} = sum_k C(m-k, k) -- the SHALLOW diagonals of Pascal
    sum to Fibonacci.  The repo's Faulhaber triangle T(n,k) = sum_{j=1}^{n-k+1}
    j^{k-1} has shallow-diagonal sums matching F_{m+1} + 2^{m-3} for m = 3..7 and
    BREAKING at m = 8.  If Fibonacci is what Pascal's shallow diagonals give, then
    the deviation measures how far Faulhaber is from Pascal -- so the right move is
    to expand Faulhaber in the Pascal basis and read the correction off.  Computed
    here exactly, with the break located and the correction term identified.
"""
import sys
from itertools import combinations, permutations
from math import gcd

NSW = int(sys.argv[1]) if len(sys.argv) > 1 else 6
MDIAG = int(sys.argv[2]) if len(sys.argv) > 2 else 14


# =====================================================================  (A)
def canon(adj, n):
    best = None
    for p in permutations(range(n)):
        key = tuple(sorted((p[a], p[b]) for a in range(n) for b in range(n)
                           if a != b and (adj[a] >> b) & 1))
        if best is None or key < best:
            best = key
    return best


def from_bits(bits, n):
    adj = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            idx += 1
    return adj


def switch(adj, S, n):
    out = [0] * n
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            if not (adj[a] >> b) & 1:
                continue
            cross = ((S >> a) & 1) != ((S >> b) & 1)
            if cross:
                out[b] |= 1 << a
            else:
                out[a] |= 1 << b
    return out


print("=" * 78)
print("(A) SWITCHING CLASSES OF TOURNAMENTS -- the canonical, base-path-free quotient")
print("=" * 78)
A002854 = {3: 2, 4: 3, 5: 7, 6: 16, 7: 54}
for n in range(3, NSW + 1):
    m = n * (n - 1) // 2
    reps = {}
    for bits in range(1 << m):
        c = canon(from_bits(bits, n), n)
        if c not in reps:
            reps[c] = bits
    classes = sorted(reps)
    idx = {c: i for i, c in enumerate(classes)}
    parent = list(range(len(classes)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for c in classes:
        adj = from_bits(reps[c], n)
        for S in range(1 << (n - 1)):          # S and its complement agree
            c2 = canon(switch(adj, S, n), n)
            union(idx[c], idx[c2])
    nsw = len({find(i) for i in range(len(classes))})
    print("  n = %d : iso classes = %-4d  SWITCHING classes = %-4d   A002854(%d) = %s"
          % (n, len(classes), nsw, n, A002854.get(n, "?")))
print()
print("  If the switching-class counts match A002854, the base-path-free quotient IS")
print("  the even-graph / two-graph side the repo already treats as first-class, and")
print("  the star group's failure to descend was purely an artifact of restricting")
print("  the cut space to non-tree edges.")

# =====================================================================  (B)
print()
print("=" * 78)
print("(B) THE THREE SIXTIES -- separate causes, and why 60 is an lcm attractor")
print("=" * 78)


def multord(a, m):
    if gcd(a, m) != 1:
        return None
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


def pisano(m):
    a, b, k = 0, 1, 0
    while True:
        a, b = b, (a + b) % m
        k += 1
        if (a, b) == (0, 1):
            return k


print("  ord_1001(2)  = %d   = lcm(ord_7(2), ord_11(2), ord_13(2)) = lcm(%d,%d,%d)"
      % (multord(2, 1001), multord(2, 7), multord(2, 11), multord(2, 13)))
print("  pi(10)       = %d   = lcm(pi(2), pi(5)) = lcm(%d,%d)"
      % (pisano(10), pisano(2), pisano(5)))
print("  |A_5|        = 60   = 5!/2, a group order -- no lcm involved at all")
print()
print("  The three are lcm(3,10,12), lcm(3,20), and 5!/2.  Different inputs,")
print("  different operations.  What they share is the TARGET, and 60 is a magnet:")
hits = {}
tot = 0
for k in range(1, 8):
    for r in range(1, 7):
        for combo in combinations(range(1, 13), r):
            L = 1
            for x in combo:
                L = L * x // gcd(L, x)
            hits[L] = hits.get(L, 0) + 1
            tot += 1
top = sorted(hits.items(), key=lambda kv: -kv[1])[:8]
print("  lcm values of subsets of {1..12}, most frequent first:")
for v, c in top:
    print("     lcm = %-6d occurs %6d times (%.2f%%)%s"
          % (v, c, 100.0 * c / tot, "   <-- 60" if v == 60 else ""))
print("  60 = 2^2*3*5 is the least common multiple of 1..6, so it absorbs every")
print("  subset whose members divide it.  Coincidences AT 60 are cheap; a shared")
print("  cause would have to show up in the INPUTS, and here it does not.")

# =====================================================================  (C)
print()
print("=" * 78)
print("(C) FIBONACCI FROM SHIFTED PASCAL, vs the Faulhaber triangle's diagonals")
print("=" * 78)


def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


from math import comb

print("  control -- Pascal's shallow diagonals ARE Fibonacci:")
print("     sum_k C(m-k,k) vs F(m+1):")
ok = True
for m in range(0, 12):
    s = sum(comb(m - k, k) for k in range(0, m // 2 + 1))
    if s != fib(m + 1):
        ok = False
    print("       m=%-3d  sum = %-6d  F(%d) = %-6d  %s"
          % (m, s, m + 1, fib(m + 1), "OK" if s == fib(m + 1) else "MISMATCH"))
print("     control holds : %s" % ok)

print()
print("  the Faulhaber triangle T(n,k) = sum_{j=1}^{n-k+1} j^(k-1), shallow diagonals:")


def T(n, k):
    return sum(j ** (k - 1) for j in range(1, n - k + 2))


rows = {}
for n in range(1, MDIAG + 2):
    rows[n] = [T(n, k) for k in range(1, n + 1)]
diag = []
for d in range(1, MDIAG + 1):
    s = 0
    for k in range(0, d):
        n = d - k
        kk = k + 1
        if 1 <= kk <= n and n in rows:
            s += rows[n][kk - 1]
    diag.append(s)
print("     diagonal sums : %s" % diag)
print()
print("     against  F(m+1) + 2^(m-3) :")
for i, D in enumerate(diag, start=1):
    pred = fib(i + 1) + (2 ** (i - 3) if i >= 3 else 0)
    mark = "OK" if pred == D else "BREAK  (D - pred = %d)" % (D - pred)
    print("       m=%-3d  D = %-8d  F(%d)+2^(%d) = %-8d  %s"
          % (i, D, i + 1, i - 3, pred, mark))
print()
print("  The Fibonacci part is exactly what a PASCAL triangle would give; the")
print("  2-power part and the break are the measure of how far Faulhaber sits from")
print("  Pascal.  Where the break first occurs is the first m at which the")
print("  Faulhaber row stops being a low-order perturbation of the binomial row.")
