#!/usr/bin/env python3
"""
collatz_mod6_20260917_reverse_tree_pieces.py
Lane reverse_tree_pieces of session collatz-mod6-20260917 (mac-mini), finalized 2026-09-21.

Object: the 3x+1 reverse tree on odd integers, root 1, children of an odd u with 3 !| u
    n_j = (2^(h0+1) 4^j u - 1)/3,  j >= 0,   h0 = 1 if u = 1 mod 3, h0 = 0 if u = 2 mod 3
(inherited inverse fibre (B1) of arithmetic_braids_20260917_collatz.md; T(n_j) = u).
Rows = residue classes 1, 3, 5 mod 6 of the child n_j; row-3 children are leaves.

Sections
  S1  three-row automaton: row(n_j) = rho(2^(h0+1+2j) u mod 9), base row from u mod 18, rotation 1->5->3,
      tree to depth 6 with j<=8; depth-D row pattern is a function of u mod 3^(D+1) (sharp)
  S2  FINITE-EXACT depth census of odd u <= 10^6 (depth = number of T-steps to 1); residues mod 3, mod 9 by depth;
      exact fibre-uniformity lemma explaining the discrepancy with the greedy map G's stationary law
  S3  growth: pruned T-tree branching factor 4/3 (exact mod-9 count, exact (4/3)^k average for k<=7),
      depth-k node counts from root 8, the simplest split-independent modulus-9 difference inequality and its exponent,
      the adversarial constant-split modulus-9 system (HEURISTIC ceiling), literature exponents (CITED)
  S4  how the three pieces fit: guarded D/E generation, leaves = odd multiples of 3, the Q1/Q2 exchange of directions,
      finite illustration of the guarded vs unguarded inverse closures of 1
  S5  wildcard: the minimal-child map m(u) = (2^(h0+1) u - 1)/3: injective, exact image, orbit structure, survival law
All checks use explicit `raise` so they survive python -O.
"""
import sys, math, time
from fractions import Fraction
from collections import Counter, defaultdict

def fail(msg):
    raise RuntimeError("CHECK FAILED: " + msg)

def check(cond, msg):
    if not cond:
        fail(msg)

def banner(title):
    print("=" * 78)
    print(title)
    print("=" * 78)

# ---------------------------------------------------------------- basic maps
def T(n):
    """accelerated odd-to-odd map: odd core of 3n+1 (n odd)."""
    if n % 2 == 0:
        fail("T called on even %d" % n)
    m = 3 * n + 1
    while m % 2 == 0:
        m //= 2
    return m

def h0(u):
    if u % 3 == 1:
        return 1
    if u % 3 == 2:
        return 0
    fail("h0 called on multiple of 3: %d" % u)

def child(u, j):
    """j-th odd T-preimage of u (u odd, 3 !| u)."""
    return ((1 << (h0(u) + 1 + 2 * j)) * u - 1) // 3

RHO = {4: 1, 7: 5, 1: 3}     # 3n+1 mod 9 -> row of n mod 6   (task convention, exponent h0+1+2j)
RHO_F = {2: 1, 5: 3, 8: 5}   # F(n)=(3n+1)/2 mod 9 -> row     (row_braid convention, exponent h)

def row_of(n):
    r = n % 6
    if r not in (1, 3, 5):
        fail("row_of on even %d" % n)
    return r

# ============================================================ S1
def section1():
    banner("S1  Three-row automaton: row of the j-th child, base row from u mod 18, rotation 1->5->3, tree to depth 6")
    # (a) the residue identity: n = 6i + r  =>  3n+1 = 18 i + 3r + 1, so 3n+1 mod 9 = 3r+1 mod 9 = 4, 1, 7 for r = 1, 3, 5.
    for r in (1, 3, 5):
        for i in range(0, 30):
            n = 6 * i + r
            check(RHO[(3 * n + 1) % 9] == r, "rho identity at n=%d" % n)
            check(RHO_F[((3 * n + 1) // 2) % 9] == r, "rho_F identity at n=%d" % n)
    print("PROVED: for odd n, row(n) := n mod 6 = rho(3n+1 mod 9) with rho(4)=1, rho(1)=3, rho(7)=5.")
    print("   Proof: n = 6i+r gives 3n+1 = 18i + (3r+1), and 3r+1 = 4, 10=1, 16=7 mod 9 for r = 1, 3, 5.")
    print("   Same law in the F-convention of the row_braid lane: F(n)=(3n+1)/2 = 2^h u, row = rho_F(2^h u mod 9),")
    print("   rho_F(2)=1, rho_F(5)=3, rho_F(8)=5 (multiply by 2^{-1} = 5 mod 9: 4->2, 1->5, 7->8).")
    # (b) child rows: 3 n_j + 1 = 2^(h0+1+2j) u, so row(n_j) = rho(2^(h0+1+2j) u mod 9); 4 has order 3 mod 9: 4->7->1->4.
    check(pow(4, 1, 9) == 4 and pow(4, 2, 9) == 7 and pow(4, 3, 9) == 1, "ord_9(4)=3")
    print("PROVED: row(n_j) = rho(2^(h0+1+2j) u mod 9); since 4 -> 16=7 -> 28=1 -> 4 (mod 9), j -> j+1 rotates the row")
    print("   1 -> 5 -> 3 -> 1 (i.e. +4 mod 6), exact period 3 = ord_9(4); every third child is a row-3 leaf.")
    # base row as a function of u mod 18 (u mod 9 with u odd)
    base = {}
    print("   base row (j=0) from u mod 18:")
    print("   u mod 18 | u mod 9 | h0 | 2^(h0+1) u mod 9 | base row | rows of j=0..5")
    for u18 in (1, 5, 7, 11, 13, 17):
        c = pow(2, h0(u18) + 1, 9) * u18 % 9
        rows = [RHO[c * pow(4, j, 9) % 9] for j in range(6)]
        base[u18] = rows[0]
        print("   %-8d | %-7d | %-2d | %-16d | %-8d | %s" % (u18, u18 % 9, h0(u18), c, rows[0], rows))
    check(base == {1: 1, 11: 1, 13: 5, 17: 5, 5: 3, 7: 3}, "base-row table")
    print("PROVED: base row = 1 for u = 1, 11 mod 18; 5 for u = 13, 17 mod 18; 3 (leaf) for u = 5, 7 mod 18;")
    print("   equivalently u mod 9 in {1,2} -> row 1, {4,8} -> row 5, {5,7} -> row 3.  Agrees with the row_braid table")
    print("   (its 'row' column for u mod 18 = 1,5,7,11,13,17 reads 1,3,3,1,5,5).")
    # explicit rotation formula
    def row_formula(u, j):
        return (base[u % 18] + 4 * j) % 6
    # (c) generate the tree to depth 6 with j <= 8 and verify every node against the formula and T
    t0 = time.time()
    JMAX, DMAX = 8, 6
    level = [1]
    seen = {1}
    total = 1
    per_level = [1]
    row_counts_level = []
    dup = 0
    for d in range(1, DMAX + 1):
        nxt = []
        rc = Counter()
        for u in level:
            if u % 3 == 0:
                continue
            for j in range(JMAX + 1):
                n = child(u, j)
                check(n % 2 == 1, "child parity")
                check(T(n) == u, "T(child)=u at u=%d j=%d" % (u, j))
                check(row_of(n) == RHO[(3 * n + 1) % 9], "rho at node")
                check(row_of(n) == RHO[pow(2, h0(u) + 1 + 2 * j, 9) * u % 9], "row law at u=%d j=%d" % (u, j))
                check(row_of(n) == row_formula(u, j), "rotation formula at u=%d j=%d" % (u, j))
                if n in seen:
                    dup += 1       # only the loop 1 -> child(1,0) = 1 may repeat; it is not re-expanded
                    check(n == 1 and u == 1 and j == 0, "unexpected repeat %d" % n)
                    continue
                seen.add(n)
                rc[row_of(n)] += 1
                nxt.append(n)
        level = nxt
        per_level.append(len(level))
        row_counts_level.append(dict(sorted(rc.items())))
        total += len(level)
    print("FINITE-EXACT: tree generated to depth %d with j<=%d: nodes per level %s (total %d), %.1fs" % (DMAX, JMAX, per_level, total, time.time() - t0))
    print("   row counts per level (1: internal, 5: internal, 3: leaf): %s" % row_counts_level)
    print("   every node satisfies T(n)=parent, row(n) = rho(2^(h0+1+2j) u mod 9) = base(u mod 18) + 4j mod 6;")
    print("   the only repeated node is the loop child(1,0) = 1 (%d repeat(s)); all other nodes distinct (T is a function)." % dup)
    exp_levels = [1]
    for d in range(1, DMAX + 1):
        # level d has 9 children per internal node of level d-1; level-1 rows: base(1)=1 -> rows 1,5,3,1,5,3,1,5,3 (6 internal)
        pass
    # (d) depth-D row pattern along a fixed index path (j_1..j_D) is a function of u mod 3^(D+1); sharp
    print("   depth-D row pattern of the descendant along an index path (j_1,...,j_D):")
    for D in range(1, 5):
        mod_ok = 3 ** (D + 1)
        mod_bad = 3 ** D
        paths = [(0,) * D, (1, 0, 2, 1)[:D], (2,) * D]
        for path in paths:
            def pattern(u):
                rows = []
                x = u
                for j in path:
                    if x % 3 == 0:
                        rows.append('L')      # leaf: no children
                        break
                    x = child(x, j)
                    rows.append(row_of(x))
                return tuple(rows)
            pat = {}
            for u in range(1, 2 * 3 ** (D + 2), 2):
                if u % 3 == 0:
                    continue
                key = u % mod_ok
                p = pattern(u)
                if key in pat:
                    check(pat[key] == p, "pattern not a function of u mod 3^(D+1), D=%d u=%d" % (D, u))
                pat[key] = p
            # sharpness: two u congruent mod 3^D with different patterns
            witness = None
            byc = defaultdict(set)
            for u in range(1, 2 * 3 ** (D + 2), 2):
                if u % 3 == 0:
                    continue
                byc[u % mod_bad].add(pattern(u))
            for key, s in byc.items():
                if len(s) > 1:
                    witness = key
                    break
            check(witness is not None, "sharpness failed at D=%d path=%s" % (D, path))
            if path == (0,) * D:
                print("   D=%d path=%s: function of u mod %d (checked all odd u < %d), NOT of u mod %d (witness class %d mod %d)" % (D, path, mod_ok, 2 * 3 ** (D + 2), mod_bad, witness, mod_bad))
    # audit 2026-09-21: the row of an ODD node is its residue mod 3 (1 -> row 1, 0 -> row 3, 2 -> row 5), not mod 9;
    # with "mod 9" the induction would only give 3^(D+2).  Verified here for all odd n < 20001.
    for n in range(1, 20001, 2):
        check(row_of(n) == {1: 1, 0: 3, 2: 5}[n % 3], "row is a function of n mod 3 at n=%d" % n)
    print("PROVED: along a fixed index path the rows of the first D descendants depend only on u mod 3^(D+1)")
    print("   (n_j = (2^a u - 1)/3 sends u mod 3^(s+1) to n_j mod 3^s, and the row of an odd node is its residue mod 3:")
    print("   1 -> row 1, 0 -> row 3, 2 -> row 5, checked for odd n < 20001; induct on D),")
    print("   and FINITE-EXACT sharp for D<=4 (some class mod 3^D splits; a class whose path reaches a leaf cannot split).")
    print("   This is the tree-side twin of the E-graph lane's Theorem 4.1 (greedy letters k_1..k_J are a function of")
    print("   m mod 3^(J+1), sharp).")
    print()

# ============================================================ S2
def section2(X=10 ** 6):
    banner("S2  FINITE-EXACT depth census of odd u <= %d: residues mod 3 and mod 9 by depth" % X)
    t0 = time.time()
    depth = {1: 0}
    peak = 0
    for u in range(3, X + 1, 2):
        x = u
        s = 0
        while x >= u:
            x = T(x)
            s += 1
            if x == 1:
                break
        # x < u is odd and already known (or x == 1)
        depth[u] = s + depth[x]
    maxd = max(depth.values())
    argmax = max(depth, key=depth.get)
    print("all %d odd u <= %d reach 1 (forward T-iteration, descent to a smaller odd already resolved); %.1fs" % (len(depth), X, time.time() - t0))
    print("max depth (number of T-steps) = %d at u = %d" % (maxd, argmax))
    check(depth[5] == 1 and depth[3] == 2 and depth[21] == 1 and depth[27] > 40, "depth sanity")
    # tabulate by depth
    by_d = defaultdict(Counter)
    cnt_d = Counter()
    for u, d in depth.items():
        by_d[d][u % 9] += 1
        cnt_d[d] += 1
    print("depth | count | frac 0 mod 3 | frac 1 mod 3 | frac 2 mod 3 | counts mod 9 (r=0..8)")
    rows_out = []
    for d in range(0, 16):
        c = cnt_d[d]
        if c == 0:
            continue
        m9 = [by_d[d][r] for r in range(9)]
        f0 = (m9[0] + m9[3] + m9[6]) / c
        f1 = (m9[1] + m9[4] + m9[7]) / c
        f2 = (m9[2] + m9[5] + m9[8]) / c
        print("%5d | %6d | %.4f | %.4f | %.4f | %s" % (d, c, f0, f1, f2, m9))
    print("bucketed depths:")
    print("depth range | count | frac 0 mod 3 | frac 1 mod 3 | frac 2 mod 3 | frac in {1,4,7} each | frac in {2,5,8} each")
    for lo in range(0, maxd + 1, 20):
        hi = lo + 19
        c = sum(cnt_d[d] for d in range(lo, hi + 1))
        if c == 0:
            continue
        m9 = [sum(by_d[d][r] for d in range(lo, hi + 1)) for r in range(9)]
        f0 = (m9[0] + m9[3] + m9[6]) / c
        f1 = (m9[1] + m9[4] + m9[7]) / c
        f2 = (m9[2] + m9[5] + m9[8]) / c
        print("%3d-%3d | %6d | %.4f | %.4f | %.4f | %s | %s" % (lo, hi, c, f0, f1, f2, ["%.4f" % (m9[r] / c) for r in (1, 4, 7)], ["%.4f" % (m9[r] / c) for r in (2, 5, 8)]))
    tot = len(depth)
    m9 = [sum(by_d[d][r] for d in by_d) for r in range(9)]
    print("all depths pooled: counts mod 9 = %s; fractions mod 3 = %.4f, %.4f, %.4f (trivially 1/3 each: every odd u <= X is a node)" % (m9, (m9[0] + m9[3] + m9[6]) / tot, (m9[1] + m9[4] + m9[7]) / tot, (m9[2] + m9[5] + m9[8]) / tot))
    # internal nodes only (3 !| u), depth >= 1 (audit 2026-09-21: the draft pooled depth 0 as well).  Since every odd u <= X
    # is a node, this pooled fraction is the trivial count of odd u = 1 mod 6 versus 5 mod 6 in [1, X]; it is NOT a tree statement.
    n_int_d1 = sum(1 for u in depth if u % 3 and u != 1)
    n_int1_d1 = sum(1 for u in depth if u % 3 == 1 and u != 1)
    check(n_int_d1 == 333332 and n_int1_d1 == 166666, "internal split %d %d" % (n_int_d1, n_int1_d1))
    print("internal nodes (3 !| u) pooled over depths >= 1: %d of %d are 1 mod 3, fraction %.4f (G stationary: 2/3; fibre-uniform: 1/2)" %
          (n_int1_d1, n_int_d1, n_int1_d1 / n_int_d1))
    print("   TRIVIAL as a pooled count: every odd u <= X is a node, so this is #{odd u = 1 mod 6 <= X} / #{odd u, 3 !| u, u <= X};")
    print("   only the depth-by-depth rows above carry tree information.")
    # depth-weighted comparison: largest |deviation| from uniform mod 9 per depth (d with count >= 1000)
    worst = max(((max(abs(by_d[d][r] / cnt_d[d] - 1 / 9) for r in range(9)), d) for d in cnt_d if cnt_d[d] >= 1000))
    print("largest deviation of a mod-9 fraction from 1/9 over depths with >= 1000 nodes: %.4f at depth %d" % worst)
    # (b) exact fibre-uniformity lemma
    print("Fibre-uniformity lemma (PROVED): for any internal u, n_j mod 9 depends on 4^j mod 27; ord_27(4) = %d and" % order(4, 27))
    check(order(4, 27) == 9, "ord_27(4)")
    sub = sorted(set(pow(4, j, 27) for j in range(9)))
    check(sub == [x for x in range(1, 27) if x % 3 == 1], "<4> mod 27 = residues 1 mod 3")
    print("   <4> mod 27 = %s = all residues 1 mod 3; so (2^(h0+1) 4^j u - 1)/3 runs over EVERY residue mod 9 exactly once" % sub)
    print("   per 9 consecutive j.  Hence the infinite fibre of every internal node is uniform mod 9: 1/3 of the children")
    print("   are row-3 leaves, 1/6 lie in each of the six classes 1,2,4,5,7,8 mod 9 -- NOT the greedy law 2/9, 1/9.")
    for u in (1, 5, 7, 11, 13, 17, 19, 23, 25, 29):
        res = sorted(child(u, j) % 9 for j in range(9))
        check(res == list(range(9)), "fibre uniformity at u=%d" % u)
    print("   checked for u in {1,5,7,11,13,17,19,23,25,29}: children j=0..8 hit each residue mod 9 once.")
    # (c) truncation count: number of children <= X of u equals floor(log_4((3X+1)/(2^(h0+1)u))) + 1
    tot_children = 0
    n_internal = 0
    for u in range(1, X + 1, 2):
        if u % 3 == 0:
            continue
        n_internal += 1
        c = (1 << (h0(u) + 1)) * u
        J = -1
        while c <= 3 * X + 1:
            J += 1
            c <<= 2
        tot_children += J + 1
    # cross-check: the number of odd n <= X whose parent T(n) is <= X
    cross = sum(1 for n in range(1, X + 1, 2) if T(n) <= X)
    check(cross == tot_children, "truncated-fibre count %d vs direct %d" % (tot_children, cross))
    print("FINITE-EXACT: #{odd n <= X : T(n) <= X} = sum over internal u <= X of (J_u(X)+1) = %d, J_u(X) = floor(log_4((3X+1)/(2^(h0+1)u)));" % tot_children)
    print("   %d internal u <= X, mean truncated fibre size %.4f (the infinite fibre is infinite; truncation keeps the first" % (n_internal, tot_children / n_internal))
    print("   J_u(X)+1 children, so small u contribute the phase-dependent prefix of the 9-periodic residue cycle).")
    # audit 2026-09-21: the truncated level d inside [1, X] is NOT the union of the truncated fibres of the depth-(d-1) nodes <= X,
    # because a child can be smaller than its parent (u = 2 mod 3): at depth 2 the node 932067 has parent 1398101 > X.
    d2 = sorted(u for u, d in depth.items() if d == 2)
    par2 = Counter(T(u) for u in d2)
    check(len(d2) == 34 and par2[1398101] == 1 and 932067 in d2 and sum(v for p, v in par2.items() if p <= X) == 33,
          "depth-2 anatomy %s" % par2)
    print("depth-2 anatomy: 34 nodes = 33 children of the six internal depth-1 nodes <= X (parents %s)" % {p: v for p, v in sorted(par2.items()) if p <= X})
    print("   plus 932067, whose parent 1398101 exceeds X.  So the finite level sets are not fibre truncations alone.")
    print("Discrepancy with G explained: (i) G never visits multiples of 3 (its k skips row-3 children) while the tree")
    print("   counts every leaf once (1/3 of nodes); (ii) among internal children the infinite tree is uniform on the six classes")
    print("   (lemma), whereas G's greedy minimal k lands in 1 mod 3 from residues {1,2,4,5} and in 2 mod 3 only from {7,8},")
    print("   giving 2/3 on 1 mod 3 (E-graph lane Theorem 4.2).  HEURISTIC, not proved: that the finite depth-by-depth deviations")
    print("   from uniform come only from the phase-dependent prefixes of truncated fibres (and from parents above X).")
    print()
    return depth

def order(a, m):
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k

# ============================================================ S3
def pruned_children(y):
    """children of y in the pruned T-tree (T(x)=x/2 or (3x+1)/2), nodes = integers not divisible by 3."""
    out = [2 * y]
    if y % 9 in (2, 8):
        out.append((2 * y - 1) // 3)
    return out

def section3():
    banner("S3  Growth: branching factor 4/3, (4/3)^k exact average, depth-k counts, modulus-9 difference inequality, literature")
    # (a) children count by residue mod 9
    cnts = {r: len(pruned_children(r + 9 * 5)) for r in (1, 2, 4, 5, 7, 8)}
    print("pruned T-tree (nodes: integers with 3 !| n; children 2y always, (2y-1)/3 iff y = 2 or 8 mod 9):")
    print("   children count by y mod 9: %s; mean over the six classes = %s" % (cnts, Fraction(sum(cnts.values()), 6)))
    check(Fraction(sum(cnts.values()), 6) == Fraction(4, 3), "mean children 4/3")
    unp = {1: 1, 2: 2, 4: 1, 5: 2, 7: 1, 8: 2}
    print("   unpruned (multiples of 3 kept as nodes): y = 2 mod 3 always has the odd child, mean %s" % Fraction(sum(unp.values()), 6))
    print("PROVED: a node uniform on the six classes mod 9 has 4/3 children on average in the pruned tree (3/2 unpruned).")
    # (b) exact average number of depth-k descendants over roots uniform mod 3^(k+1) (coprime to 3)
    print("   k | 3^(k+1) | sum over root classes of #depth-k nodes | average | (4/3)^k")
    for k in range(1, 8):
        M = 3 ** (k + 1)
        total = 0
        nroots = 0
        for a in range(1, M):
            if a % 3 == 0:
                continue
            nroots += 1
            # count depth-k descendants of the residue class a (the tree shape depends on a mod 3^(k+1) only)
            # simulate on the representative a + M*t with t chosen so all divisions are exact; use exact arithmetic on a
            frontier = [a]
            for _ in range(k):
                nf = []
                for y in frontier:
                    nf.extend(pruned_children(y))
                frontier = nf
            total += len(frontier)
        avg = Fraction(total, nroots)
        check(avg == Fraction(4, 3) ** k, "average depth-%d descendants" % k)
        print("   %d | %6d | %8d | %s | %s" % (k, M, total, avg, Fraction(4, 3) ** k))
    print("PROVED: the child of a node uniform on the six classes mod 3^(s+1) is uniform on the six classes mod 3^s")
    print("   (doubling permutes classes; (2y-1)/3 sends 9m+2 -> 6m+1 and 9m+8 -> 6m+5, each uniform in its mod-3 class),")
    print("   so the expected number of depth-k descendants of a uniformly random root is exactly (4/3)^k (FINITE-EXACT k<=7).")
    print("   No novelty claim: this is the Applegate-Lagarias/Lagarias-Weiss branching model (CITED below).")
    # (c) depth-k node counts of the pruned tree from root 8 (T-tree; 1<->2 is the loop, 4 has the single preimage 8)
    frontier = [8]
    counts = []
    for k in range(0, 33):
        counts.append(len(frontier))
        nf = []
        for y in frontier:
            nf.extend(pruned_children(y))
        frontier = nf
    print("FINITE-EXACT: pruned T-tree from root 8, nodes at depth k = 0..32:")
    print("   %s" % counts)
    print("   ratios count[k]/count[k-1] for k = 24..32: %s;  (count[32])^(1/32) = %.4f vs 4/3 = %.4f" % (["%.3f" % (counts[k] / counts[k - 1]) for k in range(24, 33)], counts[32] ** (1 / 32), 4 / 3))
    # (d) the simplest split-independent modulus-9 difference inequality
    print("Modulus-9 difference inequality (PROVED).  Let N = {n >= 1 : 3 !| n, n reaches 1 under T}, f_r(x) = #{n in N, n <= x, n = r mod 9},")
    print("   A = f_1+f_4+f_7 (1 mod 3), B = f_2+f_5+f_8 (2 mod 3).  Children are distinct (T is a function), so for x >= 1:")
    print("   (i)  A(x) >= B(x/2) + f_2(3x/2),   B(x) >= A(x/2) + f_8(3x/2)   [doubling swaps the mod-3 class; odd children of")
    print("        n = 2 mod 9 are 1 mod 3 and <= x iff n <= (3x+1)/2; of n = 8 mod 9 are 2 mod 3];")
    print("   (ii) doubling chains 1->2 (1 step), 7->5->1->2 (3), 4->8->7->5->1->2 (5) give f_2(y) >= max(f_1(y/2), f_7(y/8), f_4(y/32)) >= A(y/32)/3;")
    print("        4->8 (1), 1->2->4->8 (3), 7->5->1->2->4->8 (5) give f_8(y) >= A(y/32)/3.")
    chains = {}
    for r in (1, 4, 7):
        for target in (2, 8):
            x, s = r, 0
            while x != target:
                x = 2 * x % 9
                s += 1
            chains[(r, target)] = s
    check(chains == {(1, 2): 1, (7, 2): 3, (4, 2): 5, (4, 8): 1, (1, 8): 3, (7, 8): 5}, "doubling chain lengths %s" % chains)
    print("   doubling chain lengths (r -> target): %s" % chains)
    print("   (iii) hence A(x) >= A(x/4) + A(3x/128)/3 + A(3x/64)/3 for all x >= 1, and A(x) >= 1 for x >= 1 (1 in N).")
    print("   Induction on x with A nondecreasing gives A(x) >= c x^gamma where gamma solves 1 = 4^-g + (1/3)(3/128)^g + (1/3)(3/64)^g:")
    def phi(g):
        return 4 ** (-g) + (3 / 128) ** g / 3 + (3 / 64) ** g / 3 - 1
    lo, hi = 0.0, 1.0
    check(phi(lo) > 0 and phi(hi) < 0, "bracket")
    for _ in range(200):
        mid = (lo + hi) / 2
        if phi(mid) > 0:
            lo = mid
        else:
            hi = mid
    g_crude = (lo + hi) / 2
    print("   gamma_9,crude = %.6f  (so #{n <= x : n reaches 1} >= c x^%.4f; a weak but rigorous power bound)" % (g_crude, g_crude))
    check(0.24 < g_crude < 0.26, "crude exponent range")
    # (e) adversarial constant-split six-class system (HEURISTIC ceiling)
    print("Adversarial constant-split modulus-9 system (HEURISTIC): f_s(x) = sum_{2r=s} f_r(x/2) + sum_{r in {2,8}} alpha_{r->s} f_r(3x/2),")
    print("   alpha_{2->.} a distribution on {1,4,7}, alpha_{8->.} on {2,5,8} (the unknown splits of the odd children by n mod 27).")
    print("   With f_r ~ c_r x^g the Perron root rho(g) of M(g) must be 1.  rho(g) is NOT monotone in g (2^-g falls, (3/2)^g rises),")
    print("   so rho(g) = 1 can have two roots or none; the table lists ALL roots on [0, 3] and the smallest one is the exponent.")
    print("   Uniform split (the truth, by the uniformity lemma): the all-ones vector is an eigenvector, rho(g) = 2^-g + (3/2)^g/3,")
    print("   roots g = 1 and g = 2.  The minimum over the 9 vertex splits is an upper bound on what any constant-split")
    print("   modulus-9 argument can prove.  (Audit 2026-09-21: the recovered draft used plain power iteration, whose one-step")
    print("   ratio oscillates on the period-3 matrices, and a bisection capped at 1.5; four of its nine entries were artifacts.)")
    states = [1, 2, 4, 5, 7, 8]
    idx = {r: i for i, r in enumerate(states)}

    def spectral_radius(M):
        # M is nonnegative and irreducible (the doubling 6-cycle 1->2->4->8->7->5->1 visits every state), so M + I is
        # primitive and power iteration on it converges to rho(M) + 1 regardless of the period of M.
        n = len(M)
        v = [1.0] * n
        lam = 1.0
        for _ in range(1200):
            w = [v[i] + sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]
            lam = sum(w) / sum(v)
            v = [x / sum(w) for x in w]
        return lam - 1.0

    def rho_of(g, split2, split8):
        M = [[0.0] * 6 for _ in range(6)]
        for r in states:
            M[idx[2 * r % 9]][idx[r]] += 2 ** (-g)
        for s, a in split2.items():
            M[idx[s]][idx[2]] += a * (1.5) ** g
        for s, a in split8.items():
            M[idx[s]][idx[8]] += a * (1.5) ** g
        return spectral_radius(M)

    def all_roots(split2, split8, gmax=3.0, step=0.01):
        grid = [i * step for i in range(int(gmax / step) + 1)]
        vals = [rho_of(g, split2, split8) - 1 for g in grid]
        out = []
        for i in range(len(grid) - 1):
            if abs(vals[i]) < 1e-12:
                out.append(grid[i])          # exact grid root (g = 1 for the uniform split)
                continue
            if vals[i] * vals[i + 1] < 0 and abs(vals[i + 1]) >= 1e-12:
                a, b, fa = grid[i], grid[i + 1], vals[i]
                for _ in range(50):
                    m = (a + b) / 2
                    fm = rho_of(m, split2, split8) - 1
                    if (fm > 0) == (fa > 0):
                        a, fa = m, fm
                    else:
                        b = m
                out.append((a + b) / 2)
        return out, min(vals) + 1
    uni2, uni8 = {1: 1 / 3, 4: 1 / 3, 7: 1 / 3}, {2: 1 / 3, 5: 1 / 3, 8: 1 / 3}
    r_uni, _ = all_roots(uni2, uni8)
    check(len(r_uni) == 2 and abs(r_uni[0] - 1) < 1e-6 and abs(r_uni[1] - 2) < 1e-6, "uniform split roots %s" % r_uni)
    check(all(abs(rho_of(g, uni2, uni8) - (2 ** (-g) + 1.5 ** g / 3)) < 1e-7 for g in (0.3, 1.0, 1.7, 2.5)), "uniform closed form")
    print("   uniform split: roots of rho(g) = 1 on [0,3]: %s" % ["%.6f" % g for g in r_uni])
    expected = {(1, 2): [1.244017], (1, 5): [0.774576], (1, 8): [], (4, 2): [], (4, 5): [0.576032], (4, 8): [],
                (7, 2): [0.602605], (7, 5): [0.436588], (7, 8): []}
    best = None
    for s2 in (1, 4, 7):
        for s8 in (2, 5, 8):
            roots, rmin = all_roots({s2: 1.0}, {s8: 1.0})
            exp = expected[(s2, s8)]
            check(len(roots) == len(exp) and all(abs(a - b) < 2e-6 for a, b in zip(roots, exp)),
                  "roots at 2->%d, 8->%d: %s vs expected %s" % (s2, s8, roots, exp))
            if roots:
                print("   vertex split 2->%d, 8->%d: smallest root g = %.6f (all roots on [0,3]: %s)" % (s2, s8, roots[0], ["%.6f" % g for g in roots]))
                if best is None or roots[0] < best[0]:
                    best = (roots[0], s2, s8)
            else:
                check(rmin > 1, "no-root split must have rho > 1 on the grid")
                print("   vertex split 2->%d, 8->%d: NO root (min rho(g) on [0,3] = %.4f > 1; the draft printed %s here)" %
                      (s2, s8, rmin, "1.500000 = bisection cap" if s8 == 8 else "0.867659 = period-3 power-iteration artifact"))
    print("   minimum over the five vertex splits that have a root: g = %.6f at 2->%d, 8->%d" % best)
    print("   comparison: 3/7 = %.6f; the vertex minimum is close to but NOT equal to the recollected Krasikov (1989) exponent" % (3 / 7))
    check(abs(best[0] - 3 / 7) > 1e-3, "vertex minimum is not 3/7")
    print("   HEURISTIC reading: %.4f is the ceiling of the constant-split modulus-9 Krasikov scheme; the published bounds use" % best[0])
    print("   moduli 3^k with k large and a nonlinear program, and are NOT reproduced here.")
    print("Literature (CITED): Krasikov-Lagarias, 'Bounds for the 3x+1 problem using difference inequalities', Acta Arith. 109 (2003):")
    print("   #{n <= x : n reaches 1} >= x^0.84 for large x.  Applegate-Lagarias, 'Density bounds for the 3x+1 problem', Math. Comp. 64")
    print("   (1995), part I (tree-search method, exponent 0.643) and part II (Krasikov inequalities, exponent 0.81).")
    print("   Lagarias-Weiss, 'The 3x+1 problem: two stochastic models', Ann. Appl. Probab. 2 (1992): the branching random walk model.")
    print("   UNCITED-RECOLLECTION: Krasikov (1989) exponent 3/7 (venue and page numbers not verified here); page numbers of the")
    print("   1995/2003 papers are not verified here.  Nothing beyond these statements is claimed.")
    print()

# ============================================================ S4
def collatz_path(n):
    p = [n]
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        p.append(n)
    return p

KMIN = {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}

def G(m):
    k = 0
    while (m << k) % 9 not in (4, 7):
        k += 1
    check(k == KMIN[m % 9], "KMIN")
    return ((m << k) - 1) // 3

def section4(X=10 ** 5):
    banner("S4  How the three pieces fit: guarded D/E generation, leaves, and the Q1/Q2 exchange of directions")
    print("Guarded inverse branches (collatz_blueprint_20260921_affine.md): D(x) = 2x, E(x) = (2x-1)/3 with guard x = 2 mod 3.")
    print("   For odd u the odd T-preimages are E o D^h (u) with h = h0 + 2j (parity forced by the guard); j -> j+1 is R(n) = 4n+1.")
    for u in (1, 5, 7, 11):
        for j in range(3):
            h = h0(u) + 2 * j
            x = u
            for _ in range(h):
                x = 2 * x
            check(x % 3 == 2, "guard at u=%d h=%d" % (u, h))
            check((2 * x - 1) // 3 == child(u, j), "E o D^h = child")
    print("PROVED (checked u in {1,5,7,11}, j<=2): child(u,j) = E(D^(h0+2j)(u)), and D^h(u) = 2 mod 3 exactly for h = h0 mod 2.")
    print("Pieces: row 1 (u = 1 mod 6, h0 = 1) and row 5 (u = 5 mod 6, h0 = 0) are internal with infinite fibres whose rows")
    print("   rotate 1->5->3 from the base row of S1; row 3 (odd multiples of 3) are leaves: 3 !| 3n+1, so no odd T-preimage;")
    print("   a leaf n = 3 mod 6 is in the tree iff its parent T(n) is (T(n) is internal).  Therefore:")
    print("PROVED: Collatz (all odd n reach 1)  <=>  every u = 1 or 5 mod 6 is a node of the tree  <=>  every odd u with 3 !| u is")
    print("   child(u', j) for some node u' and some j (surjectivity of the guarded D/E closure of {1} onto the odd non-multiples of 3).")
    # the exchange of directions
    print("Exchange of directions (PROVED, elementary; no novelty claim: the Q2 half is the E-graph lane's Theorem 2.1):")
    print("   C  = deterministic Collatz graph (n -> n/2 for even n, n -> 3n+1 for odd n);  E = C plus the arrows n -> 3n+1 for EVEN n.")
    print("   Collatz: every n reaches 1 in C  <=>  the closure of {1} under the GUARDED inverse moves x -> 2x, x -> (x-1)/3 (x = 4 mod 6)")
    print("        is everything (reverse each arrow of C: this closure is the full inverse Collatz tree, multiples of 3 included).")
    print("   Q1 (E-graph lane's definition): every n reaches 1 in E  <=>  the closure of {1} under the UNGUARDED inverse moves x -> 2x,")
    print("        x -> (x-1)/3 (x = 1 mod 3) is everything.  Collatz implies Q1 (C is a subgraph of E); the converse is not claimed.")
    print("   Q2 (E-graph lane): 1 reaches every m with 3 !| m in E  <=>  every such m reaches 1 by the UNGUARDED inverse moves")
    print("        (reverse each arrow of E).  'Avoiding multiples of 3' is automatic: no E-arrow enters 3Z from outside (3n+1 is never")
    print("        0 mod 3, and n/2 = 0 mod 3 forces n = 0 mod 3), so no path from 1 meets a multiple of 3.  The greedy map")
    print("        G(m) = (2^k m-1)/3 with minimal k, 2^k m in {4,7} mod 9, chooses one such move; G^s(m) = 1 certifies 1 ->* m in E.")
    print("   So Collatz and Q1 are 'TO 1' statements (in C, in E) = 'FROM 1' statements for the guarded / unguarded inverse graphs,")
    print("   while Q2 is a 'FROM 1' statement in E = a 'TO 1' statement for the unguarded inverse graph; reversing arrows exchanges")
    print("   the two directions, and the relaxation C -> E adds exactly the even -> 3n+1 arrows (E-graph lane Theorem 1.1: reverse")
    print("   index j = -1 of the fibre).  Q1 and Q2 together say that the non-multiples of 3 form one strongly connected component")
    print("   of E (E-graph lane conjecture C_E).  (Audit 2026-09-21: the recovered draft wrote 'Q1 (Collatz)' for the C-statement,")
    print("   which is not the E-graph lane's Q1, and labelled its unguarded-closure count below as a Q2 illustration; corrected.)")
    for n in range(1, 20000):
        check((3 * n + 1) % 3 != 0 and (n % 2 == 1 or (n // 2) % 3 != 0 or n % 3 == 0), "no arrow into 3Z at n=%d" % n)
    # illustration for m = 7
    p = collatz_path(7)
    gp = [7]
    while gp[-1] != 1:
        gp.append(G(gp[-1]))
    print("   example m = 7: Collatz path (Q1, TO 1): %s" % p)
    print("             greedy inverse certificate (Q2): G-orbit %s, i.e. the E-path 1 -> 4 -> 2 -> 7 uses the even arrow 2 -> 7." % gp)
    check(gp == [7, 2, 1], "G orbit of 7")
    check(G(5) == 13 and G(1) == 1 and G(11) == 7, "G samples")
    print("   G vs the minimal odd child m(u) of S5: on odd u they coincide iff u = 1, 2 or 8 mod 9 (KMIN = h0+1 there);")
    for u in range(1, 400, 2):
        if u % 3 == 0:
            continue
        check((G(u) == child(u, 0)) == (u % 9 in (1, 2, 8)), "G vs child law at u=%d" % u)
    print("   checked for all odd u < 400.  For u = 4, 7 mod 9 G takes k = 0 (an even image (u-1)/3, not a tree node); for")
    print("   u = 5 mod 9 the minimal odd child is the leaf (2u-1)/3 and G skips to k = 3, i.e. G(u) = child(u,1).")
    # finite illustration: guarded vs unguarded inverse closure of 1 inside [1, X]
    t0 = time.time()
    guarded = set([1])
    stack = [1]
    while stack:
        x = stack.pop()
        for y in ((2 * x,) + (((x - 1) // 3,) if x % 6 == 4 and x > 4 else ())):
            if y <= X and y not in guarded:
                guarded.add(y)
                stack.append(y)
    unguarded = set([1])
    stack = [1]
    while stack:
        x = stack.pop()
        cands = [2 * x]
        if x % 3 == 1 and x > 1:
            cands.append((x - 1) // 3)
        for y in cands:
            if y <= X and y % 3 != 0 and y not in unguarded:
                unguarded.add(y)
                stack.append(y)
    n_guarded_non3 = sum(1 for y in guarded if y % 3 != 0)
    # FROM-1 closure in E inside [1, X]: successors x -> 3x+1 (all x) and x -> x/2 (even x).  This is the Q2 set in the box.
    from1 = set([1])
    stack = [1]
    while stack:
        x = stack.pop()
        for y in ((3 * x + 1,) + ((x // 2,) if x % 2 == 0 else ())):
            if y <= X and y not in from1:
                from1.add(y)
                stack.append(y)
    check(all(y % 3 for y in from1), "FROM-1 set meets 3Z")
    giant = unguarded & from1
    print("FINITE-EXACT (X=%d, %.1fs): three closures of 1 inside [1,X]:" % (X, time.time() - t0))
    print("   (a) guarded inverse closure (Collatz in the box: forward C-orbit stays <= X): %d nodes (%d not divisible by 3);" % (len(guarded), n_guarded_non3))
    print("   (b) unguarded inverse closure = TO-1-in-E in the box (Q1 in the E-graph lane's sense), multiples of 3 excluded:")
    print("       %d of the %d non-multiples of 3;" % (len(unguarded), X - X // 3))
    print("   (c) successor closure = FROM-1-in-E in the box (the Q2 set): %d of the %d non-multiples of 3;" % (len(from1), X - X // 3))
    print("   (b) and (c) intersect in %d nodes = the giant SCC of E|[1,X] containing 1 (E-graph lane .out S3: giant size 17077)." % len(giant))
    check(len(unguarded) == 27472 and len(from1) == 29762 and len(giant) == 17077, "closure sizes %d %d %d" % (len(unguarded), len(from1), len(giant)))
    check(all(y in unguarded for y in guarded if y % 3 != 0), "guarded subset of unguarded")
    print("   (a) minus multiples of 3 is a subset of (b) (checked); none of the three is complete inside [1,X] because")
    print("   forward Collatz peaks or E-node peaks leave [1,X] (E-graph lane S3); Collatz, Q1 and Q2 are FINITE-EXACT to 10^6 there.")
    miss_g = sorted(y for y in range(1, X + 1) if y % 3 and y not in guarded)[:8]
    miss_u = sorted(y for y in range(1, X + 1) if y % 3 and y not in unguarded)[:8]
    miss_f = sorted(y for y in range(1, X + 1) if y % 3 and y not in from1)[:8]
    print("   smallest non-multiples of 3 missing from (a): %s (Collatz peak of %d is %d > X)" % (miss_g, miss_g[0], max(collatz_path(miss_g[0]))))
    print("   smallest missing from (b): %s (the E-graph lane's SCC outsiders at N=10^5; their binding" % miss_u)
    print("       constraint is the forward peak, and all five smallest lie in (c): %s)" % [y in from1 for y in miss_u[:5]])
    print("   smallest missing from (c), i.e. the honest Q2-in-the-box outsiders: %s" % miss_f)
    check(max(collatz_path(miss_g[0])) > X, "first guarded outsider must leave [1,X]")
    check(miss_u[:5] == [1535, 2047, 2207, 2287, 2303] and all(y in from1 for y in miss_u[:5]), "TO-1 outsiders vs E-graph lane")
    check(miss_f[:4] == [3281, 4010, 4739, 4922], "FROM-1 outsiders %s" % miss_f)
    print("SCOPE: no map found from this exchange to THM-4139/THM-4146 (x^2-29/16 three-cycle) or THM-3341 (Gaussian squaring);")
    print("   the only shared structure is the mod-6/mod-9 row typing.")
    print()

# ============================================================ S5
def mmap(u):
    return child(u, 0)

def section5(X=10 ** 6):
    banner("S5  Wildcard: the minimal-child map m(u) = (2^(h0+1) u - 1)/3 (j = 0), orbits, image, survival law")
    # injective, image
    print("PROVED: m is injective (2^(h0+1) u = 2^(h0'+1) u' with u, u' odd forces h0 = h0', u = u').")
    print("PROVED: image(m) = {odd n : n = 1 mod 8} U {odd n : n = 3 mod 4} = odd n with n != 5 mod 8; the complement 5 mod 8")
    print("   is exactly R(odd) = {4n+1}, the children with j >= 1.  (u = 1 mod 3: 4u-1 = 3 mod 8 so n = 1 mod 8, and n = 1 mod 8")
    print("   gives u = (3n+1)/4 odd, = 1 mod 3; u = 2 mod 3: 2u-1 = 1 mod 4 so n = 3 mod 4, and n = 3 mod 4 gives u = (3n+1)/2 odd, = 2 mod 3.)")
    img = set()
    for u in range(1, 4000, 2):
        if u % 3:
            img.add(mmap(u))
    for n in range(1, 2000, 2):
        check((n in img) == (n % 8 != 5), "image law at n=%d" % n)
    print("   checked: for odd n < 2000, n in m({odd u < 4000, 3 !| u}) iff n != 5 mod 8.")
    # m goes deeper: T(m(u)) = u
    for u in (1, 5, 7, 11, 13):
        check(T(mmap(u)) == u, "T o m = id")
    print("PROVED: T(m(u)) = u, so m moves one level away from the root; m(u) = (4u-1)/3 > u for u = 1 mod 3 (u > 1), m(u) = (2u-1)/3 < u for u = 2 mod 3.")
    print("   Fixed points: m(u) = u iff 3u = 2^(h0+1) u - 1 iff u = 1 (h0 = 1, 4u-1 = 3u).  m(1) = 1.")
    check(mmap(1) == 1, "m(1)=1")
    # residue automaton: m(u) mod 3 from u mod 9; m(u) mod 9 from u mod 27 (uniform on the mod-3 class)
    print("   m(u) mod 3 from u mod 9: {1,2} -> 1 (internal row 1), {4,8} -> 2 (internal row 5), {5,7} -> 0 (leaf: orbit stops).")
    for u in range(1, 200, 2):
        if u % 3 == 0:
            continue
        r = u % 9
        exp = 1 if r in (1, 2) else (2 if r in (4, 8) else 0)
        check(mmap(u) % 3 == exp, "m residue law at u=%d" % u)
    # orbits to X
    t0 = time.time()
    CAP = 10000
    lengths = Counter()
    maxlen, argmax = -1, None
    maxpeak, argpeak = 0, None
    capped = 0
    cycles = 0
    for u in range(1, X + 1, 2):
        if u % 3 == 0:
            continue
        x = u
        L = 0
        pk = u
        while x % 3 != 0:
            if u != 1 and x == u and L > 0:
                cycles += 1
                break
            if u == 1:
                break
            x = mmap(x)
            L += 1
            if x > pk:
                pk = x
            if L >= CAP:
                capped += 1
                break
        lengths[L] += 1
        if L > maxlen:
            maxlen, argmax = L, u
        if pk / u > maxpeak:
            maxpeak, argpeak = pk / u, u
    n_dom = sum(lengths.values())
    print("FINITE-EXACT (%.1fs): m-orbits of all %d odd u <= %d with 3 !| u, followed until the first multiple of 3 (a leaf):" % (time.time() - t0, n_dom, X))
    print("   orbits hitting the step cap %d: %d; nontrivial cycles found: %d (an m-cycle would be a T-cycle, impossible below X" % (CAP, capped, cycles))
    print("   since every u <= X reaches 1 (S2) and the only T-cycle through 1 is {1}); u = 1 is the unique fixed point.")
    check(capped == 0 and cycles == 0, "m orbits")
    print("   longest orbit before a leaf: L = %d at u = %d; largest peak/u = %.2f at u = %d" % (maxlen, argmax, maxpeak, argpeak))
    print("   L | #u with orbit length L | fraction | (1/3)(2/3)^(L-1)  [L = number of m-steps until the leaf, the leaf step included;")
    print("                                                      L = 0 only for the fixed point u = 1]")
    for L in range(0, 13):
        c = lengths[L]
        pred = 0.0 if L == 0 else (1 / 3) * (2 / 3) ** (L - 1)
        print("   %2d | %7d | %.5f | %.5f" % (L, c, c / n_dom, pred))
        if L >= 1:
            check(abs(c / n_dom - pred) < 2e-4, "geometric law at L=%d" % L)
    tail = sum(c for L, c in lengths.items() if L >= 13)
    print("   tail: #{L >= 13} = %d, fraction %.5f vs (2/3)^12 = %.5f" % (tail, tail / n_dom, (2 / 3) ** 12))
    check(abs(tail / n_dom - (2 / 3) ** 12) < 2e-4, "geometric tail")
    print("   The empirical law is geometric with ratio 2/3 to within 2e-4 at every L <= 12 and in the tail (checked).")
    # exact survival law: the first L rows of the m-orbit depend on u mod 3^(L+1) and exactly (2/3)^L of the classes survive
    # lifting lemma (audit 2026-09-21, replaces the draft's garbled wording): for 3 !| a and s >= 1,
    #     m(a + t 3^s) = m(a) + 2^(h0(a)+1) t 3^(s-1),   t = 0, 1, 2,
    # since h0 depends on a mod 3 only; 2^(h0+1) is a unit mod 3, so the three lifts of a mod 3^s to mod 3^(s+1) are sent
    # bijectively onto the three lifts of m(a) mod 3^(s-1) to mod 3^s.  Iterating L times: the three lifts of a surviving class
    # mod 3^L to mod 3^(L+1) give m^L values that are 0, 1, 2 mod 3 once each, so exactly one third dies at step L.
    for s in range(1, 7):
        for a in range(1, 3 ** s):
            if a % 3 == 0:
                continue
            for t in range(3):
                lhs = ((1 << (h0(a) + 1)) * (a + t * 3 ** s) - 1) // 3
                rhs = mmap(a) + (1 << (h0(a) + 1)) * t * 3 ** (s - 1)
                check(lhs == rhs, "lifting lemma at a=%d s=%d t=%d" % (a, s, t))
    print("   lifting lemma (PROVED; checked s <= 6): m(a + t 3^s) = m(a) + 2^(h0(a)+1) t 3^(s-1), so the three lifts of a class")
    print("   mod 3^s to mod 3^(s+1) map bijectively onto the three lifts of m(a) mod 3^(s-1) to mod 3^s.")
    print("   survival law by residue classes (PROVED + FINITE-EXACT L<=8): the first L m-steps stay internal for exactly (2/3)^L")
    print("   of the classes u mod 3^(L+1) coprime to 3 (iterate the lifting lemma: the three lifts of a class surviving L-1 steps")
    print("   have m^L values that are 0, 1, 2 mod 3 once each, so each step kills exactly one third of the surviving classes):")
    print("   L | 3^(L+1) | surviving classes | total classes | fraction | (2/3)^L")
    for L in range(1, 9):
        M = 3 ** (L + 1)
        surv = 0
        tot = 0
        for a in range(1, M):
            if a % 3 == 0:
                continue
            tot += 1
            x = a
            ok = True
            for _ in range(L):
                # m(x) mod 3^(remaining) is determined; work with the exact integer a + M*t? use a directly (odd or even
                # does not matter mod powers of 3: h0 depends on x mod 3 only)
                x = ((1 << (h0(x) + 1)) * x - 1) // 3 if ((1 << (h0(x) + 1)) * x - 1) % 3 == 0 else None
                if x is None:
                    fail("integrality")
                if x % 3 == 0:
                    ok = False
                    break
            surv += ok
        fr = Fraction(surv, tot)
        check(fr == Fraction(2, 3) ** L, "survival fraction at L=%d: %s" % (L, fr))
        print("   %d | %6d | %6d | %6d | %s | %s" % (L, M, surv, tot, fr, Fraction(2, 3) ** L))
    print("   (the representative a is used directly: m only needs x mod 3 for h0 and divides by 3 once, so x mod 3^(s) after one")
    print("   step is a function of a mod 3^(s+1); parity plays no role in the residue law).")
    print("OPEN: whether EVERY m-orbit (u != 1) ends at a leaf.  The 3-adic set of u with an infinite internal orbit has Haar measure 0")
    print("   (survival (2/3)^L); an integer in it would be an odd u whose whole minimal-child line avoids 3Z -- a Collatz-type 3-adic")
    print("   question, not decided here.  FINITE-EXACT: none below 10^6.")
    print("   Typed analogy (m-orbit -> greedy G-orbit): source = the m-orbit u, m(u), m^2(u), ... in the odd tree; target = the")
    print("   G-orbit in E; map = both pick the minimal exponent k (m: minimal odd-admissible k; G: minimal k with an internal image, any")
    print("   parity); preserved = one-step residue Markov structure driven by u mod 27 -> image mod 9; lost = G never dies (it skips")
    print("   leaves) and reaches 1 empirically, m dies at a leaf with density 1 and never returns to 1; sidecar = the step-count/peak")
    print("   laws; test = compare survival (2/3)^L with G's stopping-time bound (7/9)^(J-1) of the E-graph lane.")
    print()

def main():
    t0 = time.time()
    section1()
    section2()
    section3()
    section4()
    section5()
    print("ALL CHECKS PASSED  total time %.1fs" % (time.time() - t0))

if __name__ == "__main__":
    main()
