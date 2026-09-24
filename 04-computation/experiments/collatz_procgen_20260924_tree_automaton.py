#!/usr/bin/env python3
"""collatz_procgen_20260924_tree_automaton.py

Lane: inverse tree mod 192 (session collatz-procgen-20260922, 2026-09-24).

Part 1 of the owner's programme, made exact.

Shortcut map T(n) = n/2 (n even), (3n+1)/2 (n odd).  Inverse moves:
    D(x) = 2x                      (always legal)
    E(x) = (2x-1)/3                (legal iff x = 2 mod 3; then E(x) is odd)
    S(p) = 4p+1                    (sibling ladder; E D^2 = S E on legal points)
The minus sheet T_-(n) = n/2, (3n-1)/2 has E_-(x) = (2x+1)/3, legal iff x = 1 mod 3.

Sections (printed):
  A1  well-definedness of the moves on types mod 3*2^k (k = 1..10)
  A2  the owner's automaton A_6 on Z/192: SCCs, condensation, branching, spectral radius
  A3  finite-level shadows: D-attractor, E-self-loop at -1, S-cycle {21,85,149}, root 2-cycle
  A4  what a residue mod 192 decides (parity window, mod-3 types, T^6 mod 3^(1+a), predecessors)
  A5  order relations inside a 6-step window: exact gate bound and census on both sheets
  A6  the owner's dictionary: siblings, 4x recursion, (N-1)/3 vs (4N-1)/3, phantom sibling = E-graph arrow
  A7  equivalence theorem checks (tree of {1,2} = basin of {1,2}; unicyclic functional graph)
  A8  the negation isomorphism nu(r) = -r between the plus and minus automata, k = 1..12

Every check raises on failure.  Runtime about 2 s; peak memory about 300 MB.
"""
import sys
import time
from fractions import Fraction
from math import comb, log, gcd

import numpy as np

T0 = time.time()


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)


def T(n, b=1):
    return n // 2 if n % 2 == 0 else (3 * n + b) // 2


def D(x):
    return 2 * x


def E_legal(x, b=1):
    # E_b(x) = (2x - b)/3 legal iff 3 | (2x - b)
    return (2 * x - b) % 3 == 0


def E(x, b=1):
    assert E_legal(x, b)
    return (2 * x - b) // 3


def S(p):
    return 4 * p + 1


def v2(n):
    n = abs(n)
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


def mem(tag):
    import resource
    print(f"[mem] {tag}: max RSS so far {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB", file=sys.stderr)


def hdr(s):
    mem("before " + s[:3])
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# ---------------------------------------------------------------------------
hdr("A1  Well-definedness of the moves on types mod 3*2^k")
# ---------------------------------------------------------------------------
print("Claim F: T(n) mod 3*2^(k-1) is a function of n mod 3*2^k.")
print("Claim D: 2n mod 3*2^(k+1) is a function of n mod 3*2^k (injective, image = even classes).")
print("Claim E: for n = 2 mod 3, E(n) mod 2^(k+1) is a function of n mod 3*2^k; the map is a")
print("         bijection onto the odd classes mod 2^(k+1); E(n) mod 3 is NOT a function of")
print("         n mod 3*2^k: the three lifts mod 9*2^k give E(n) = 0,1,2 mod 3 once each.")
print("Claim E': E(n) mod 3*2^j is a function of n mod 9*2^(j-1) (j >= 1), and not of n mod 3*2^k for any k.")
for b in (1, -1):
    for k in range(1, 11):
        M = 3 * 2 ** k
        LIFTS = 6
        # F
        for r in range(M):
            vals = {T(r + M * t, b) % (M // 2) for t in range(LIFTS)}
            check(len(vals) == 1, f"F well-defined b={b} k={k} r={r}")
        # D
        img = set()
        for r in range(M):
            vals = {(2 * (r + M * t)) % (2 * M) for t in range(LIFTS)}
            check(len(vals) == 1, f"D well-defined k={k}")
            img |= vals
        check(img == {x for x in range(2 * M) if x % 2 == 0}, "D image = even classes")
        # E
        guard = 2 if b == 1 else 1
        img = []
        for r in range(M):
            if r % 3 != guard:
                continue
            vals = {E(r + M * t, b) % (2 ** (k + 1)) for t in range(LIFTS)}
            check(len(vals) == 1, f"E mod 2^(k+1) well-defined b={b} k={k}")
            img.append(vals.pop())
            mod3 = [E(r + M * t, b) % 3 for t in range(3)]
            check(sorted(mod3) == [0, 1, 2], f"E mod 3 free over the three lifts b={b} k={k} r={r}")
        check(sorted(img) == list(range(1, 2 ** (k + 1), 2)), f"E bijective onto odd classes b={b} k={k}")
    # E' : E(n) mod 3*2^j from n mod 9*2^(j-1)
    for j in range(1, 9):
        M9 = 9 * 2 ** (j - 1)
        for r in range(M9):
            if r % 3 != (2 if b == 1 else 1):
                continue
            vals = {E(r + M9 * t, b) % (3 * 2 ** j) for t in range(6)}
            check(len(vals) == 1, f"E' b={b} j={j}")
print("  verified for both sheets, k = 1..10 (Claims F, D, E) and j = 1..8 (Claim E'),")
print("  six lifts per class.  In particular E(n) mod 192 needs n mod 288 = 9*32, D(n) mod 192")
print("  needs n mod 96, and T(n) mod 96 needs n mod 192.")

# ---------------------------------------------------------------------------
hdr("A2  The owner's inverse-tree type automaton A_k at modulus 3*2^k; k = 6 (192) in detail")
# ---------------------------------------------------------------------------
print("States Z/3*2^k.  Edges:  r -D-> 2r mod 3*2^k  (every r);")
print("                         r -E-> s  for the three s with s = (2r-1)/3 mod 2^k   (r = 2 mod 3).")
print("The E-edge is nondeterministic exactly in the mod-3 digit (A1, Claim E); each of the three")
print("targets carries Haar weight 1/3.  This is a sound over-approximation: m in T^-1(n) implies")
print("type(m) in delta(type(n)).")


def build_inverse_automaton(k, b=1):
    M = 3 * 2 ** k
    P = 2 ** k
    guard = 2 if b == 1 else 1
    edges = {r: [] for r in range(M)}
    for r in range(M):
        edges[r].append(("D", (2 * r) % M))
        if r % 3 == guard:
            e2 = E(r, b) % P  # well defined mod 2^(k+1), so mod 2^k
            for s3 in range(3):
                # s = e2 mod 2^k and s = s3 mod 3 (CRT)
                s = next(s for s in range(e2, M, P) if s % 3 == s3)
                edges[r].append(("E", s))
    return edges


def tarjan(nodes, succ):
    index = {}
    low = {}
    onstack = set()
    stack = []
    sccs = []
    counter = [0]
    sys.setrecursionlimit(100000)

    def strong(v):
        index[v] = low[v] = counter[0]
        counter[0] += 1
        stack.append(v)
        onstack.add(v)
        for w in succ[v]:
            if w not in index:
                strong(w)
                low[v] = min(low[v], low[w])
            elif w in onstack:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = stack.pop()
                onstack.discard(w)
                comp.append(w)
                if w == v:
                    break
            sccs.append(sorted(comp))

    for v in nodes:
        if v not in index:
            strong(v)
    return sccs


def scc_report(k, b=1, verbose=False):
    M = 3 * 2 ** k
    edges = build_inverse_automaton(k, b)
    succ = {r: sorted({s for (_, s) in edges[r]}) for r in range(M)}
    sccs = tarjan(range(M), succ)
    sizes = sorted((len(c) for c in sccs), reverse=True)
    big = max(sccs, key=len)
    comp_of = {}
    for i, c in enumerate(sccs):
        for v in c:
            comp_of[v] = i
    nontriv = [c for c in sccs if len(c) > 1 or c[0] in succ[c[0]]]
    return edges, succ, sccs, sizes, big, comp_of, nontriv


print()
print(" k  | states | SCC sizes (largest first, then counts)            | big SCC = non-multiples of 3? | nontrivial SCCs")
for k in range(1, 11):
    edges, succ, sccs, sizes, big, comp_of, nontriv = scc_report(k)
    M = 3 * 2 ** k
    nonmult = [r for r in range(M) if r % 3 != 0]
    same = (big == nonmult)
    check(same, f"big SCC is exactly the non-multiples of 3 at k={k}")
    check(sizes[0] == 2 ** (k + 1) and all(s == 1 for s in sizes[1:]) and len(sizes) == 1 + 2 ** k,
          f"SCC sizes at k={k}")
    check(len(nontriv) == 2 and [0] in nontriv, f"nontrivial SCCs = big and {{0}} at k={k}")
    print(f" {k:2d} | {M:6d} | [{sizes[0]}] + {len(sizes)-1} singletons                              | {same}                          | big, {{0}} (D-self-loop)")

k = 6
M = 192
edges, succ, sccs, sizes, big, comp_of, nontriv = scc_report(6)
# reachability from the root types {1,2}
reach = set()
frontier = [1, 2]
while frontier:
    v = frontier.pop()
    if v in reach:
        continue
    reach.add(v)
    frontier.extend(succ[v])
check(reach == set(range(M)), "every type mod 192 reachable from the root types {1,2}")
# mult-of-3 region: closed, D only, a DAG into 0
mult3 = [r for r in range(M) if r % 3 == 0]
for r in mult3:
    check(all(s % 3 == 0 for s in succ[r]) and len(succ[r]) == 1, "multiples of 3 are closed under D only")
depth_to0 = {}
for r in mult3:
    x, d = r, 0
    while x != 0:
        x = (2 * x) % M
        d += 1
        check(d <= 6, "D nilpotent on the 2-part")
    depth_to0[r] = d
hist = {}
for r in mult3:
    hist[depth_to0[r]] = hist.get(depth_to0[r], 0) + 1
print()
print("A_6 (modulus 192):")
print(f"  SCCs: one of size {len(big)} = the 128 non-multiples of 3 (2-part free, mod-3 part in {{1,2}}),")
print(f"        and 64 singletons = the multiples of 3; only {{0}} carries a loop (D(0) = 0).")
print(f"  Condensation: big SCC --E(choice 0 mod 3)--> odd multiples of 3 --D-chain--> 0 (sink).")
print(f"  Steps to reach 0 inside the multiples of 3 (D is nilpotent on the 2-part mod 64): {dict(sorted(hist.items()))}")
print(f"  Every one of the 192 types is reachable from the root types {{1, 2}}.")
outdeg = {}
for r in range(M):
    outdeg[len(edges[r])] = outdeg.get(len(edges[r]), 0) + 1
print(f"  Out-degree census (with multiplicity): {dict(sorted(outdeg.items()))}   (64 branch states r = 2 mod 3 have D + 3 E-choices)")

# weighted matrix: D weight 1, each E-choice weight 1/3
W = np.zeros((M, M))
for r in range(M):
    for (lab, s) in edges[r]:
        W[r, s] += 1.0 if lab == "D" else 1.0 / 3.0
rows = W.sum(axis=1)
mean_children = rows.mean()
check(abs(mean_children - 4 / 3) < 1e-12, "mean weighted out-degree 4/3")
ev = np.linalg.eigvals(W)
rho = max(abs(ev))
check(abs(rho - 4 / 3) < 1e-9, "spectral radius of the weighted inverse automaton = 4/3")
core = [r for r in range(M) if r % 3 != 0]
Wc = W[np.ix_(core, core)]
evc = np.linalg.eigvals(Wc)
top = sorted(evc, key=lambda z: -abs(z))[:4]
# right eigenvector for 4/3 on the core: depends on r mod 3 only, ratio v2/v1 = 4/3
vals, vecs = np.linalg.eig(Wc)
i43 = int(np.argmin(abs(vals - 4 / 3)))
v = np.real(vecs[:, i43])
v = v / v[core.index(1)]
v1 = v[[core.index(r) for r in core if r % 3 == 1]]
v2_ = v[[core.index(r) for r in core if r % 3 == 2]]
check(np.allclose(v1, 1.0) and np.allclose(v2_, 4 / 3), "Perron vector (1, 4/3) on classes 1, 2 mod 3")
print(f"  Weighted matrix (D weight 1, each E-choice 1/3): mean row sum = {mean_children:.6f} = 4/3;")
print(f"  spectral radius = {rho:.12f} = 4/3; on the core the Perron vector is 1 on r = 1 mod 3 and 4/3 on r = 2 mod 3")
print(f"  (solves 3L^2 - L - 4 = 0, roots 4/3 and -1); leading core eigenvalues: " +
      ", ".join(f"{z.real:+.6f}{z.imag:+.6f}i" for z in top))
print("  Pruned branching (children that are not multiples of 3): from r = 1 mod 3: 1 (the D-child);")
print("  from r = 2 mod 3: 1 + 2/3; average over the core 4/3.  Probability that an E-child is a leaf: 1/3.")

# forward automaton on 192: r -> {T(r), T(r+192)} mod 192
fsucc = {r: sorted({T(r) % M, T(r + M) % M}) for r in range(M)}
fsccs = tarjan(range(M), fsucc)
fbig = max(fsccs, key=len)
check(fbig == core, "forward automaton: big SCC = non-multiples of 3")
check(len(fsccs) == 1 + 64, "forward automaton: 64 transient multiples of 3")
for r in mult3:
    for s in fsucc[r]:
        pass
# multiples of 3 cannot be re-entered from the core
for r in core:
    check(all(s % 3 != 0 for s in fsucc[r]), "forward: core never enters 3Z")
print()
print("  Forward automaton on Z/192 (r -> T(r) mod 192, two choices of the top bit): one SCC of the")
print("  128 non-multiples of 3 and 64 singleton SCCs (the multiples of 3; only 0 carries a loop, T(0) = 0).")
print("  The core never enters 3Z: after the first odd step an orbit never meets 3Z again.  So 3Z is a")
print("  SOURCE region forward and a SINK region backward: the leaf cones of the inverse tree.")

# ---------------------------------------------------------------------------
hdr("A3  Finite-level shadows of the generators and of the root at modulus 192")
# ---------------------------------------------------------------------------
# D dynamics mod 192
Dmap = {r: (2 * r) % M for r in range(M)}
periodic_D = set()
for r in range(M):
    x = r
    for _ in range(20):
        x = Dmap[x]
    y = x
    cyc = [y]
    while True:
        y = Dmap[y]
        if y == x:
            break
        cyc.append(y)
    periodic_D.add(tuple(sorted(cyc)))
check(periodic_D == {(0,), (64, 128)}, "D-periodic types mod 192")
# E self-loop at -1 = 191
check(191 % 3 == 2 and E(191) % 128 == 127, "E(-1) = -1: E-self-loop at 191")
check(("E", 191) in edges[191], "191 -E-> 191 in A_6")
# S dynamics mod 192
Smap = {r: (4 * r + 1) % M for r in range(M)}
periodic_S = set()
for r in range(M):
    x = r
    for _ in range(20):
        x = Smap[x]
    cyc = [x]
    y = Smap[x]
    while y != x:
        cyc.append(y)
        y = Smap[y]
    periodic_S.add(tuple(sorted(cyc)))
check(periodic_S == {(21, 85, 149)}, "S-periodic types mod 192 = {21,85,149}")
check(all((3 * r + 1) % 64 == 0 for r in (21, 85, 149)), "21,85,149 = -1/3 mod 64")
# root cycle
check(("D", 2) in edges[1] and ("E", 1) in edges[2], "root 2-cycle 1 -D-> 2 -E-> 1")
trunk = [(4 ** i - 1) // 3 for i in range(1, 30)]
check(all(trunk[i - 1] == E(D(1) if i == 1 else 2 ** (2 * i - 1)) for i in range(1, 30)), "trunk = E D^(2i-1)(1)")
check(all(trunk[i] == S(trunk[i - 1]) for i in range(1, 29)), "trunk = S^(i-1)(1)")
tr192 = [t % M for t in trunk]
print("  D: fixed type 0 and the 2-cycle {64, 128} (the 2-part is nilpotent, the mod-3 digit swaps 1<->2).")
print("     Real meaning: 0 is the fixed point of D; deep D-chains of a multiple of 3 end in type 0.")
print("  E: the type 191 = -1 mod 192 has an E-self-loop, E(-1) = -1 exactly (E's fixed point is -1).")
print("  S: every sibling ladder p, 4p+1, 16p+5, ... ends in the 3-cycle {21, 85, 149} mod 192:")
print("     all three residues mod 3 of -1/3 mod 64 (3*21+1 = 64).  S^j(p) -> -1/3 in Z_2, and -1/3 = E(0).")
print("  Root: 1 -D-> 2 -E-> 1 is a 2-cycle of A_6 (the cycle {1,2} of T).")
print(f"  Trunk (4^i-1)/3 = E D^(2i-1)(1) = S^(i-1)(1) mod 192: {tr192[:9]} ... (enters {{21,85,149}} at i = 3).")

# ---------------------------------------------------------------------------
hdr("A4  What a residue mod 192 decides")
# ---------------------------------------------------------------------------


def orbit_word(n, L, b=1):
    xs = [n]
    w = []
    for _ in range(L):
        w.append(xs[-1] % 2)
        xs.append(T(xs[-1], b))
    return xs, w


LIFT = 3 ** 8  # lifts r + 192 t, t < 3^8, cover every class mod 3^9 * 64
dec_rows = []
for b in (1, -1):
    for r in range(M):
        base_xs, base_w = orbit_word(r if r > 0 else M, 6, b)
        a = sum(base_w)
        types = [x % 3 for x in base_xs]
        t6 = base_xs[6] % 3 ** (a + 1)
        for t in range(1, 60):
            n = r + M * t
            xs, w = orbit_word(n, 6, b)
            check(w == base_w, "parity window x_0..x_5 decided by n mod 64")
            check([x % 3 for x in xs] == types, "types mod 3 of T^0..T^6 decided by n mod 192")
            check(xs[6] % 3 ** (a + 1) == t6, "T^6(n) mod 3^(1+a) decided by n mod 192")
        # T^6 mod 3^(a+2) is NOT decided (show for one lift pair)
        if b == 1:
            vals = {orbit_word(r + M * t, 6, b)[0][6] % 3 ** (a + 2) for t in range(1, 12)}
            check(len(vals) == 3, "T^6(n) mod 3^(a+2) takes 3 values over lifts")
        dec_rows.append((b, r, a))
print("  For every r mod 192 and both sheets (59 lifts each):")
print("   * n mod 64 decides the parity window (x_0, ..., x_5) of the T-orbit (Terras; a bijection Z/64 -> {0,1}^6);")
print("   * n mod 192 decides the mod-3 types of T^0(n), ..., T^6(n): after the first odd step the type is")
print("     2 at an odd-step image and alternates 1,2 along halvings; before it, it is 2^(-j) n mod 3;")
print("   * n mod 192 decides T^6(n) mod 3^(1+a), a = #odd steps in the window, and not mod 3^(2+a);")
print("   * backward: n mod 192 decides D(n) mod 384 and, for n = 2 mod 3, E(n) mod 128, but not E(n) mod 3.")
print("  Read on a tree node m = w(root): m mod 64 is the LAST six moves of its path (E = odd, D = even),")
print("  m mod 3 is its own branching type.  So the owner's type is: one branching trit + a six-move window.")
# distribution of a over classes
ahist = {}
for (b, r, a) in dec_rows:
    if b == 1:
        ahist[a] = ahist.get(a, 0) + 1
check(all(ahist[a] == 3 * comb(6, a) for a in range(7)), "odd-step counts binomial")
print(f"  #classes mod 192 by odd steps a in the window: {dict(sorted(ahist.items()))} = 3*C(6,a).")

# ---------------------------------------------------------------------------
hdr("A5  Order relations inside the 6-step window: exact gate bound and census, both sheets")
# ---------------------------------------------------------------------------
print("For a window x_i -> x_j (j-i = m <= 6) with a odd steps: 2^m x_j = 3^a x_i + b*c, c >= 0 (c = 0 iff a = 0).")
print("Word prediction: x_j > x_i iff 3^a > 2^m.  A crossing (actual order != prediction) needs")
print("x_i < gate = c/|2^m - 3^a| on the side where b*c and 2^m - 3^a have the same sign.")


def words(m):
    for bits in range(2 ** m):
        yield [(bits >> t) & 1 for t in range(m)]


def carry(w, b):
    # affine map of the word: x -> (3^a x + b*c)/2^m ; return a, c
    A, C, den = 1, 0, 1  # x_t = (A x + C)/den
    for bit in w:
        if bit:
            A, C = 3 * A, 3 * C + b * den
            den *= 2
        else:
            den *= 2
    # x_m = (A x + C)/den with den = 2^m
    return A, C


gates = {}
for b in (1, -1):
    G = Fraction(0)
    argG = None
    for m in range(1, 7):
        for w in words(m):
            A, C = carry(w, b)
            den = 2 ** m
            if A == den:
                continue
            # x_j - x_i = ((A - den) x + C)/den ; prediction sign(A - den)
            # crossing iff sign((A-den)x + C) != sign(A-den) for x>0 : needs C*(A-den) < 0 and x < |C|/|A-den|
            if C * (A - den) < 0:
                g = Fraction(abs(C), abs(A - den))
                if g > G:
                    G, argG = g, (m, tuple(w))
    gates[b] = (G, argG)
    print(f"  sheet b={b:+d}: largest gate over all words of length <= 6 = {G} = {float(G):.4f} at (m, word) = {argG}")
check(gates[1][0] == Fraction(76, 5) and gates[-1][0] == Fraction(260, 17), "6-window gate bounds 76/5 and 260/17")
# census on n <= 10^6: non-generic windows
N = 10 ** 6
for b in (1, -1):
    bad_starts = set()
    bad_windows = 0
    n = np.arange(1, N + 1, dtype=np.int64)
    xs = [n]
    for _ in range(6):
        x = xs[-1]
        xs.append(np.where(x % 2 == 0, x // 2, (3 * x + b) // 2))
    par = [(x % 2).astype(np.int64) for x in xs[:6]]
    for i in range(6):
        a = np.zeros(N, dtype=np.int64)
        for j in range(i + 1, 7):
            a = a + par[j - 1]
            m = j - i
            pred_up = (3 ** a) > (2 ** m)
            actual_up = xs[j] > xs[i]
            eq = xs[j] == xs[i]
            badmask = (pred_up != actual_up) | eq
            if badmask.any():
                idx = np.nonzero(badmask)[0]
                bad_windows += len(idx)
                bad_starts.update((n[idx]).tolist())
                check(np.all(xs[i][idx] <= float(gates[b][0])), "every crossing starts below the gate bound")
    bs = sorted(bad_starts)
    print(f"  sheet b={b:+d}: n <= 10^6 with a non-generic comparison inside the window T^0..T^6: {len(bs)} starts,"
          f" {bad_windows} windows (equalities at cycle points included); starts = {bs[:40]}")
    check(max(bs) <= 64 * float(gates[b][0]), "all non-generic starts are small")
print("  PROVED: every value in the window is >= n/64 (each step at least halves), so for n > 64*76/5 (plus)")
print("  resp. n > 64*260/17 (minus), i.e. for every n > 979, the full order pattern of (T^0 n, ..., T^6 n)")
print("  is the word's generic pattern: a function of n mod 64, the same on both sheets after r -> -r.")
print("  The census shows the true exceptions are the listed starts (all <= 32 plus, <= 80 minus).")
print("  (Consistent with the order-laws note: short-lag genericity; the first sheet-separating window")
print("  has 12 odd steps, far beyond the 6 steps a residue mod 192 sees.)")

# ---------------------------------------------------------------------------
hdr("A6  The owner's dictionary: siblings, the 4x recursion and the hidden choice")
# ---------------------------------------------------------------------------
# preimages
for n in range(1, 20001):
    pre = sorted({m for m in (2 * n, (2 * n - 1) // 3) if m >= 1 and T(m) == n})
    exp = [2 * n] + ([E(n)] if n % 3 == 2 else [])
    check(sorted(exp) == pre, f"T^-1({n})")
    if n % 3 == 2:
        check(E(n) % 2 == 1 and E(n) < n, "E(n) odd and smaller")
# brute-force: no other preimage below 3n
for n in range(1, 1001):
    pre = [m for m in range(1, 3 * n + 3) if T(m) == n]
    exp = sorted([2 * n] + ([E(n)] if n % 3 == 2 else []))
    check(pre == exp, "brute-force preimages")
print("  T^-1(n) = {2n} U {(2n-1)/3 : n = 2 mod 3} (brute force n <= 1000, identity n <= 20000).")
# E D^2 = S E
for x in range(-3000, 3001):
    if x % 3 == 2:
        check(E(D(D(x))) == S(E(x)), "E D^2 = S E")
        check(E_legal(4 * x), "legality preserved by D^2")
print("  E D^2 = S E with S(p) = 4p+1 on every legal x in [-3000, 3000] (identity: (8x-1)/3 = 4(2x-1)/3 + 1).")


def U(p):  # Syracuse map on odd p
    q = 3 * p + 1
    while q % 2 == 0:
        q //= 2
    return q


for p in range(1, 20001, 2):
    check(U(S(p)) == U(p), "siblings share the Syracuse image")
    check(v2(3 * S(p) + 1) == v2(3 * p + 1) + 2, "v2 ladder")
    check(3 * S(S(p)) + 1 == 16 * (3 * p + 1), "3 S^j(p) + 1 = 4^j (3p+1)")
print("  U(4p+1) = U(p), v_2(3S(p)+1) = v_2(3p+1)+2 and 3S^j(p)+1 = 4^j(3p+1) (odd p < 20000):")
print("  the siblings converge 2-adically to -1/3 = E(0), whose T-image is 0 = D's fixed point.")
# (N-1)/3 and (4N-1)/3
cnt_true = cnt_phantom = 0
for N_ in range(1, 30001):
    if N_ % 3 != 1:
        continue
    lo = (N_ - 1) // 3
    hi = (4 * N_ - 1) // 3
    check(hi == S(lo), "(4N-1)/3 = S((N-1)/3)")
    if N_ % 2 == 0:
        check(lo % 2 == 1 and T(lo) == N_ // 2 and T(hi) == 2 * N_, "genuine siblings when N even")
        check(E(N_ // 2) == lo and E(2 * N_) == hi, "both are E-images")
        cnt_true += 1
    else:
        check(lo % 2 == 0 and 3 * lo + 1 == N_, "phantom sibling = even x with 3x+1 = N (E-graph arrow)")
        check(T(hi) == 2 * N_, "(4N-1)/3 still a genuine T-predecessor of 2N")
        cnt_phantom += 1
print(f"  For N = 1 mod 3 (N <= 30000): (4N-1)/3 = 4*(N-1)/3 + 1, consecutive rungs of one sibling ladder.")
print(f"   N even ({cnt_true} cases): (N-1)/3 = E(N/2) is odd, a genuine T-predecessor of N/2; the 'A/2 branch'.")
print(f"   N odd  ({cnt_phantom} cases): (N-1)/3 is EVEN; it is the phantom rung j = -1, and the arrow")
print("          (N-1)/3 -> 3*(N-1)/3 + 1 = N is exactly the extra arrow of the E-graph relaxation.")
print("  So the owner's 'hidden choice between both odds' is the sibling index; A/2 exists iff the lower rung")
print("  is a genuine T-odd; the relaxation E-SCC is the ladder extended one rung down.")

# ---------------------------------------------------------------------------
hdr("A7  The equivalence theorem, checked")
# ---------------------------------------------------------------------------
print("Theorem 1 (PROVED in the note): R := closure of {1,2} under D and guarded E equals the T-basin of {1,2};")
print("hence Collatz (every n >= 1 reaches 1) iff R = Z_{>0}.  Checks:")
Nf = 10 ** 6
# forward: every n <= Nf reaches {1,2}
reach1 = np.zeros(Nf + 1, dtype=bool)
reach1[1] = reach1[2] = True
for n0 in range(3, Nf + 1):
    x = n0
    while x >= n0:
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
    reach1[n0] = reach1[x]
check(reach1[1:].all(), "every n <= 10^6 reaches {1,2}")
print(f"  (i) every n <= {Nf} reaches {{1,2}} under T (descent to an already-resolved smaller value).")
# backward closure with a value cap
X = 3000
CAP = 10 ** 5
seen = {1, 2}
stack = [1, 2]
while stack:
    x = stack.pop()
    for y in ((2 * x,) + ((E(x),) if x % 3 == 2 else ())):
        if y <= CAP and y not in seen:
            seen.add(y)
            stack.append(y)
inR = {n for n in range(1, X + 1) if n in seen}
# forward: n <= X whose orbit stays <= CAP
stay = set()
for n0 in range(1, X + 1):
    x, ok = n0, True
    while x not in (1, 2):
        x = T(x)
        if x > CAP:
            ok = False
            break
    if ok:
        stay.add(n0)
check(inR == stay, "capped closure = capped basin")
print(f"  (ii) the D/E-closure of {{1,2}} inside [1, 10^5] contains exactly the n <= {X} whose T-orbit")
print(f"       stays <= 10^5 ({len(inR)} of {X}; the others peak above the cap, e.g. {sorted(set(range(1, X+1)) - inR)[:6]}).")
check(len(inR) < X, "the cap is binding for some n (non-monotone sizes in the inverse tree)")
# unicyclic: T maps R to R, the only cycle is {1,2}
check(T(1) == 2 and T(2) == 1 and E(2) == 1 and D(1) == 2, "root cycle edges")
print("  (iii) functional graph: every node has exactly one parent T(n); by (i) the only cycle meeting")
print("        [1, 10^6] is {1,2} (E(2) = 1 and D(1) = 2 are its two edges).")

# ---------------------------------------------------------------------------
hdr("A8  The negation isomorphism nu(r) = -r between the plus and minus automata")
# ---------------------------------------------------------------------------
for k in range(1, 13):
    Mk = 3 * 2 ** k
    ep = build_inverse_automaton(k, 1)
    em = build_inverse_automaton(k, -1)
    for r in range(Mk):
        mapped = sorted((lab, (-s) % Mk) for (lab, s) in ep[r])
        check(mapped == sorted(em[(-r) % Mk]), f"nu maps plus edges to minus edges k={k} r={r}")
        # forward type map
    # forward map check via integers of both signs
    for n in range(-5 * Mk, 5 * Mk):
        check(T(-n, -1) == -T(n, 1), "T_-(-n) = -T_+(n)")
print("  For k = 1..12 the bijection r -> -r mod 3*2^k maps every D- and E_+-edge of the plus automaton")
print("  onto a D- resp. E_- -edge of the minus automaton (guard 2 mod 3 <-> 1 mod 3), and T_-(-n) = -T_+(n).")
print("  Root types: {1,2} (plus) <-> {-1,-2} (the minus sheet's NEGATIVE cycle); the minus sheet's positive")
print("  fixed point 1 <-> the plus type -1, i.e. the E_+-fixed point.")

print()
print(f"[tree_automaton] ALL CHECKS PASSED   ({time.time() - T0:.1f} s)", file=sys.stderr)
print("[tree_automaton] ALL CHECKS PASSED")
