#!/usr/bin/env python3
"""collatz_mod6_20260922_w6_sign_specific_probes.py

Lane sign_specific_probes (session collatz-mod6-20260922, wave 6).

Six candidate sign-specific invariants are run against the 3n-1 sheet
(T_-(n) = (3n-1)/2^v on positive odd n; three cycles with odd skeletons
(1), (5,7), (17,25,37,55,41,61,91)) and the 3n+1 sheet, and each is sorted
honestly into: reduces to the sign law, sheet-blind, or genuinely new.

Sections (each prints its own header):
  S0  session-lead probe: square-sum graph Q_n components, degree<=1 vertices,
      Hamiltonian paths n<=32; the two pasted errors (3,7,11,17 not an AP,
      the pasted Lean adjacency is not loopless)
  S1  (i)   the root asymmetry in the summand reading, halved and unhalved,
            and the companion family (a n + b)/2^k; strictness census to 10^6
  S2  (ii)  Berggren transport: (3,+-1,k) edges with hypotenuse <= 10^6, 10^7;
            the O(log X) sheet-difference bound; where the sheets touch
  S3  (iii) E_- reachability: Q1_- and Q2_- to 10^6, one-SCC picture
  S4  (iv)  carry statistics B_L mod 4, mod 3^L, v_2 on cycles and orbits
  S5  (v)   greedy G_+-, cycle words vs greedy reverse words, inverse 2-cycles
  S6  (vi)  the word-function theorem and the order-only residue-of-minimum law
All checks use explicit `raise`; timing goes to stderr only.
"""
import math
import sys
import time
from fractions import Fraction

T0 = time.time()


def tlog(msg):
    sys.stderr.write("[%7.2fs] %s\n" % (time.time() - T0, msg))


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# ---------------------------------------------------------------------------
# shared maps
# ---------------------------------------------------------------------------
def v2(n):
    if n == 0:
        raise ValueError("v2(0)")
    return (n & -n).bit_length() - 1


def T(n, b):
    """odd-skeleton step: (3n+b)/2^k, returns (image, k)."""
    m = 3 * n + b
    k = v2(m)
    return m >> k, k


def word_of(n, b, L):
    w = []
    for _ in range(L):
        n, k = T(n, b)
        w.append(k)
    return tuple(w)


def carry(word, b):
    """2^K T^L(n) = 3^L n + B, B = b * sum 3^(L-1-i) 2^(K_i), K_0 = 0."""
    L = len(word)
    B = 0
    K = 0
    for i, k in enumerate(word):
        B += 3 ** (L - 1 - i) * 2 ** K
        K += k
    return b * B, K


def greedy_k(m, b):
    """minimal k >= 0 with (2^k m - b)/3 an integer not divisible by 3."""
    for k in range(0, 8):
        t = (m << k) - b
        if t % 3 == 0 and (t // 3) % 3 != 0:
            return k
    raise RuntimeError("no greedy k for m=%d b=%d" % (m, b))


def G(m, b):
    k = greedy_k(m, b)
    return ((m << k) - b) // 3, k


# ---------------------------------------------------------------------------
# S0  session-lead probe: square-sum graph
# ---------------------------------------------------------------------------
hdr("S0  session-lead probe: square-sum graph Q_n (x~y iff x+y square, x!=y)")


def sq_adj(n):
    sqs = set(i * i for i in range(2, int(math.isqrt(2 * n)) + 2))
    adj = {x: [] for x in range(1, n + 1)}
    for x in range(1, n + 1):
        for y in range(x + 1, n + 1):
            if x + y in sqs:
                adj[x].append(y)
                adj[y].append(x)
    return adj


def components(adj):
    seen = set()
    comps = []
    for s in adj:
        if s in seen:
            continue
        stack = [s]
        seen.add(s)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))
    return comps


def ham_path(adj, n, budget=3_000_000):
    """backtracking Hamiltonian path; returns (path or None, nodes, exhausted)."""
    low = [x for x in adj if len(adj[x]) <= 1]
    if any(len(adj[x]) == 0 for x in adj) and n > 1:
        return None, 0, False
    if len(low) >= 3:
        return None, 0, False
    starts = low if low else list(adj)
    nodes = [0]
    exhausted = [False]
    visited = set()
    path = []

    def rec(u):
        nodes[0] += 1
        if nodes[0] > budget:
            exhausted[0] = True
            return False
        path.append(u)
        visited.add(u)
        if len(path) == n:
            return True
        # unvisited neighbours ordered by remaining degree (Warnsdorff)
        cand = [w for w in adj[u] if w not in visited]
        cand.sort(key=lambda w: sum(1 for z in adj[w] if z not in visited))
        for w in cand:
            if rec(w):
                return True
            if exhausted[0]:
                return False
        path.pop()
        visited.discard(u)
        return False

    for s in starts:
        path.clear()
        visited.clear()
        if rec(s):
            return list(path), nodes[0], False
        if exhausted[0]:
            return None, nodes[0], True
    return None, nodes[0], False


print("n | #components | degree<=1 vertices | Hamiltonian path")
ham_yes = []
ham_no = []
comp_counts = {}
for n in range(1, 33):
    adj = sq_adj(n)
    comps = components(adj)
    low = [x for x in adj if len(adj[x]) <= 1]
    p, nodes, exh = ham_path(adj, n)
    if exh:
        raise RuntimeError("Hamiltonian search budget exhausted at n=%d" % n)
    comp_counts[n] = len(comps)
    if p is not None:
        for i in range(n - 1):
            s = p[i] + p[i + 1]
            check(math.isqrt(s) ** 2 == s, "bad path step")
        check(sorted(p) == list(range(1, n + 1)), "path not Hamiltonian")
        ham_yes.append(n)
    else:
        ham_no.append(n)
    print("%2d | %2d | %-14s | %s" % (n, len(comps), str(low) if low and n >= 14 else "-",
                                     "yes " + ",".join(map(str, p)) if p else "no"))
print("Hamiltonian path exists for n in", ham_yes)
print("no Hamiltonian path for n in", ham_no)
check(ham_yes == [1, 15, 16, 17, 23] + list(range(25, 33)), "Hamiltonian list vs session lead")
check(all(comp_counts[n] == 3 for n in range(4, 13)), "3 components 4<=n<=12")
check(comp_counts[13] == 2, "2 components at 13")
check(all(comp_counts[n] == 1 for n in range(14, 33)), "connected 14..32")
for n, expect in [(18, [16, 17, 18]), (19, [16, 18])] + [(n, [18]) for n in range(20, 31)] + [(31, []), (32, [])]:
    adj = sq_adj(n)
    low = [x for x in adj if len(adj[x]) <= 1]
    check(low == expect, "degree<=1 set at n=%d: %s" % (n, low))
lead15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
lead23 = [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22]
for path in (lead15, lead23):
    n = len(path)
    check(sorted(path) == list(range(1, n + 1)), "lead path not a permutation")
    for i in range(n - 1):
        s = path[i] + path[i + 1]
        check(math.isqrt(s) ** 2 == s, "lead path step not square")
    sq_used = sorted(set(path[i] + path[i + 1] for i in range(n - 1)))
    print("session-lead path n=%d verified; squares used: %s" % (n, sq_used))
print("squares used by the n=15 path:", sorted(set(lead15[i] + lead15[i + 1] for i in range(14))),
      "(the 'braid 9,16,25' of the paste is exact for this path)")
print("session lead citation (not re-fetched here): OEIS A090461 = 15,16,17,23,25,26,... (every k >= 25 per its comments)")
print("pasted 'Delta=4 ladder' 3,7,11,17: gaps", [7 - 3, 11 - 7, 17 - 11], "-> not an arithmetic progression")
print("pasted Lean Adj x y := exists k, (x+1)+(y+1)=k^2 at x=y=1 (value 2): 2+2=4=2^2 -> a loop; not loopless")
# in this note's Q_n with x!=y the only square 2x is excluded; loops would be x with 2x square: 2,8,18,32
print("values x<=32 with 2x a perfect square (loops in the pasted definition):",
      [x for x in range(1, 33) if math.isqrt(2 * x) ** 2 == 2 * x])
tlog("S0 done")

# ---------------------------------------------------------------------------
# S1  (i) the root asymmetry in the summand reading
# ---------------------------------------------------------------------------
hdr("S1  (i) root asymmetry: summand companions of the odd arrows, both sheets")
print("halved reading   (T-form, k=1 arrow): n -> n + c,  c_+ = (n+1)/2,  c_- = (n-1)/2")
print("unhalved reading (C/E-form arrow):    n -> n + c,  c_+ = 2n+1,     c_- = 2n-1")
print("strict summand arrow needs: companion c positive and c != n")
print()
print("n | halved c_+ status | halved c_- status | unhalved c_+ status | unhalved c_- status")


def status(c, n):
    if c <= 0:
        return "companion %d not positive: NOT AN ARROW" % c
    if c == n:
        return "companion %d = n: DIAGONAL (excluded)" % c
    return "strict"


for n in [1, 3, 5, 7, 9]:
    print("%d | %s | %s | %s | %s" % (n, status((n + 1) // 2, n), status((n - 1) // 2, n),
                                      status(2 * n + 1, n), status(2 * n - 1, n)))
# census to 10^6 of the four predicates
X1 = 10 ** 6
fail = {"h+": [], "h-": [], "u+": [], "u-": []}
for n in range(1, X1 + 1, 2):
    c = (n + 1) // 2
    if c <= 0 or c == n:
        fail["h+"].append(n)
    c = (n - 1) // 2
    if c <= 0 or c == n:
        fail["h-"].append(n)
    c = 2 * n + 1
    if c <= 0 or c == n:
        fail["u+"].append(n)
    c = 2 * n - 1
    if c <= 0 or c == n:
        fail["u-"].append(n)
print()
print("odd n <= %d where the arrow is not a strict summand arrow:" % X1)
for key, name in [("h+", "halved plus"), ("h-", "halved minus"), ("u+", "unhalved plus"), ("u-", "unhalved minus")]:
    print("  %-14s : %s" % (name, fail[key]))
check(fail["h+"] == [1] and fail["h-"] == [1] and fail["u+"] == [] and fail["u-"] == [1], "root census")
print("PROVED (algebra): (n+1)/2 = n iff n = 1; (n-1)/2 <= 0 iff n = 1; 2n+1 = n never (n>0); 2n-1 = n iff n = 1.")
print("So the strict-summand shadow differs between the sheets ONLY at n = 1, and in a reading-dependent way:")
print("  halved:   plus 1 -> 2 is the excluded diagonal (1,1); minus 1 -> 1 has companion 0 (not an arrow)")
print("  unhalved: plus 1 -> 4 = 1 + 3 is strict;              minus 1 -> 2 = 1 + 1 is the excluded diagonal")
print()
print("companion family: the k=1 companion map n -> (n+b)/2 is the a=1 member of (a n + b)/2^k;")
print("1 is a fixed point of (a n + b)/2^k iff a + b is a power of two:")
for a, b in [(3, 1), (3, -1), (1, 1), (1, -1), (5, 1), (5, -1), (7, 1), (7, -1)]:
    s = a + b
    fixed = s > 0 and (s & (s - 1)) == 0
    print("  (a,b)=(%d,%2d): a+b=%2d -> 1 fixed: %s" % (a, b, s, fixed))
check((1 + 1) & 1 == 0 and (1 - 1) == 0, "family root")
# the minus cycles' arrows are all strict summand arrows (k=1 steps) or descents (k>=2)
cycles_minus = [(1,), (5, 7), (17, 25, 37, 55, 41, 61, 91)]
cycle_plus = [(1,)]
print()
print("cycle arrows, halved reading (k=1 steps are summand arrows n -> n + (n+b)/2; k>=2 steps are descents):")
for b, cyc in [(1, cycle_plus[0])] + [(-1, c) for c in cycles_minus]:
    L = len(cyc)
    parts = []
    for i in range(L):
        n = cyc[i]
        m, k = T(n, b)
        check(m == cyc[(i + 1) % L], "cycle step")
        if k == 1:
            c = (n + b) // 2
            parts.append("%d->%d (k=1, companion %d, %s)" % (n, m, c, "strict" if 0 < c != n else ("diagonal" if c == n else "NOT AN ARROW")))
        else:
            parts.append("%d->%d (k=%d, %s)" % (n, m, k, "descent" if m < n else "fixed point, companion 0"))
    print("  sheet %+d cycle %s: %s" % (b, cyc, "; ".join(parts)))
tlog("S1 done")

# ---------------------------------------------------------------------------
# S2  (ii) Berggren transport: hypotenuse counts of (3,+-1,k) edges
# ---------------------------------------------------------------------------
hdr("S2  (ii) (3,+-1,k) edges x->y (3x+b = 2^k y, x,y odd, x != y) with hypotenuse (x^2+y^2)/2 <= X")


def edges_upto(X, b):
    """directed edges x->y with (x^2+y^2)/2 <= X; returns dict k -> list of (x,y)."""
    out = {}
    xmax = math.isqrt(2 * X) + 1
    for x in range(1, xmax + 1, 2):
        m = 3 * x + b
        if m <= 0:
            continue
        k = v2(m)
        y = m >> k
        if y == x:
            continue
        if (x * x + y * y) <= 2 * X:
            out.setdefault(k, []).append((x, y))
    return out


def summarize(X):
    res = {}
    for b in (1, -1):
        ed = edges_upto(X, b)
        tot = sum(len(v) for v in ed.values())
        tri = set()
        for k, lst in ed.items():
            for x, y in lst:
                tri.add((max(x, y), min(x, y)))
        res[b] = (ed, tot, len(tri))
    return res


C_EDGE = sum(math.sqrt(2.0 / (1.0 + 9.0 / 4 ** k)) / 2 ** (k + 1) for k in range(1, 60))
print("C_edge = sum_k sqrt(2/(1+9/4^k))/2^(k+1) = %.6f (sheet-blind main term)" % C_EDGE)
for X in (10 ** 6, 10 ** 7):
    res = summarize(X)
    edp, totp, trip = res[1]
    edm, totm, trim = res[-1]
    ks = sorted(set(edp) | set(edm))
    print()
    print("X = %d: directed edges plus %d (triangles %d), minus %d (triangles %d), C_edge*sqrt(X) = %.1f"
          % (X, totp, trip, totm, trim, C_EDGE * math.sqrt(X)))
    print("  k | #plus | #minus | diff")
    for k in ks:
        print("  %2d | %5d | %6d | %+d" % (k, len(edp.get(k, [])), len(edm.get(k, [])), len(edp.get(k, [])) - len(edm.get(k, []))))
    kmax = max(ks)
    print("  |plus - minus| = %d <= 2*(k_max+1) = %d (k_max = %d)" % (abs(totp - totm), 2 * (kmax + 1), kmax))
    check(abs(totp - totm) <= 2 * (kmax + 1), "difference bound")
    if X == 10 ** 6:
        check(totp == 507 and totm == 506 and trim == 505, "berggren lane 507/506/505 reproduction")
        union = set()
        for ed in (edp, edm):
            for lst in ed.values():
                for x, y in lst:
                    union.add((max(x, y), min(x, y)))
        print("  union of triangles over b=+-1: %d from %d directed edges (b=+1 and b=-1 triangle sets disjoint: %s)"
              % (len(union), totp + totm, len(union) == trip + trim))
        check(len(union) == 1012, "union 1012")
    print("  ratio plus/minus = %.6f" % (totp / totm))

# where do the sheets touch in the Berggren tree?  Generic B3 of a k=1 edge of sheet b is a
# k=2 edge (2x+b)->y of sheet -b (cited: berggren lane Thm 3.1, control S6).  Count them within X.
print()
print("sheet-mixing parent->child pairs (B3 child of a k=1 edge lands on the other sheet), hypotenuse <= X:")
for X in (10 ** 6, 10 ** 7):
    res = summarize(X)
    for b in (1, -1):
        ed = res[b][0]
        k1 = ed.get(1, [])
        n_child_in = 0
        for x, y in k1:
            u = 2 * x + b
            m = 3 * u - b
            kk = v2(m)
            check((m >> kk) == y and kk == 2, "generic B3 identity")
            if u * u + y * y <= 2 * X:
                n_child_in += 1
        print("  X=%d sheet %+d: %d k=1 edges, %d of their B3 children (k=2 edges of sheet %+d) also have hypotenuse <= X"
              % (X, b, len(k1), n_child_in, -b))
print("  so the sheets touch at EVERY k=1 edge, not only in the root cluster {(3,1),(5,1),(5,3),(13,5),(7,5),(19,7),(9,7)};")
print("  the root cluster is the unique legal component with more than three pairs (cited Thm 3.2), i.e. the only")
print("  place where mixing chains to depth > 1.")
# smallest mixing pairs listing
res6 = summarize(2000)
print("  first ten B3 mixing pairs on each sheet (x->y k=1 of sheet b) => ((2x+b)->y k=2 of sheet -b):")
for b in (1, -1):
    lst = sorted(res6[b][0].get(1, []))[:10]
    print("   b=%+d:" % b, ", ".join("%d->%d => %d->%d" % (x, y, 2 * x + b, y) for x, y in lst))
tlog("S2 done")

# ---------------------------------------------------------------------------
# S3  (iii) E_- reachability to 10^6
# ---------------------------------------------------------------------------
hdr("S3  (iii) E_- (arrows n -> n/2 for even n, n -> 3n-1 for all n): Q1_- and Q2_- to 10^6")
X3 = 10 ** 6
# T_- basin census on odd n <= X3 (memoised by first value below n)
basin = {1: 1, 5: 5, 17: 17}
cyc_of = {}
for cyc in cycles_minus:
    for c in cyc:
        cyc_of[c] = cyc[0]
maxsteps = 0
for n in range(1, X3 + 1, 2):
    if n in cyc_of:
        basin[n] = cyc_of[n]
        continue
    m = n
    steps = 0
    while m >= n and m not in cyc_of:
        m, _ = T(m, -1)
        steps += 1
        if steps > 100000:
            raise RuntimeError("step cap at n=%d" % n)
    maxsteps = max(maxsteps, steps)
    basin[n] = cyc_of[m] if m in cyc_of else basin[m]
cnt = {1: 0, 5: 0, 17: 0}
for n in range(1, X3 + 1, 2):
    cnt[basin[n]] += 1
print("T_- basin census, odd n <= %d: basin{1} %d, basin{5,7} %d, basin{17..91} %d, escapes 0, max steps to drop %d"
      % (X3, cnt[1], cnt[5], cnt[17], maxsteps))
check(sum(cnt.values()) == X3 // 2, "basin total")
# explicit E_- escape paths from the cycle minima to 1
esc5 = [5, 14, 7, 20, 10, 29, 86, 43, 128, 64, 32, 16, 8, 4, 2, 1]
esc17 = [17, 50, 25, 74, 37, 110, 55, 164, 82, 41, 122, 61, 182, 91, 272, 136, 68, 203, 608, 304, 911, 2732, 1366, 683,
         2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1]


def is_E_minus_arrow(u, w):
    return (u % 2 == 0 and w == u // 2) or w == 3 * u - 1


for path in (esc5, esc17):
    for i in range(len(path) - 1):
        check(is_E_minus_arrow(path[i], path[i + 1]), "E_- arrow %d -> %d" % (path[i], path[i + 1]))
    evens_used = [(path[i], path[i + 1]) for i in range(len(path) - 1) if path[i] % 2 == 0 and path[i + 1] == 3 * path[i] - 1]
    print("E_- path %d -> 1 with %d arrows verified; even->3n-1 arrows used: %s" % (path[0], len(path) - 1, evens_used))
    print("  path: " + " -> ".join(map(str, path)))
print("Q1_- (every n reaches 1 in E_-): FINITE-EXACT for n <= %d: each n halves to an odd n' <= n, n' reaches a cycle" % X3)
print("  under T_- (census above, T_- steps are E_- paths), and the three cycle minima reach 1 by the paths above.")
# Q2_- greedy certification to 10^6 (reproduces control S4)
starts = 0
below = 0
hit_one = 0
cycled = []
capped = 0
hist = {}
for m in range(2, X3 + 1):
    if m % 3 == 0:
        continue
    starts += 1
    x = m
    seen = set()
    st = 0
    status_ = None
    while True:
        x, _ = G(x, -1)
        st += 1
        if x < m:
            status_ = "below"
            if x == 1:
                hit_one += 1
            break
        if x in seen:
            status_ = "cycle"
            break
        seen.add(x)
        if st > 2000:
            status_ = "cap"
            break
    if status_ == "below":
        below += 1
        hist[st] = hist.get(st, 0) + 1
    elif status_ == "cycle":
        cycled.append(m)
    else:
        capped += 1
print("Q2_- greedy G_-(m) = (2^k m + 1)/3, k minimal with result an integer not 0 mod 3, starts 2 <= m <= %d, 3 !| m: %d"
      % (X3, starts))
print("  certified below m: %d; hit 1 first: %d; cycled: %s; capped: %d" % (below, hit_one, cycled, capped))
print("  steps-to-certify histogram (first six):", [(s, hist[s]) for s in sorted(hist)[:6]], "max steps", max(hist))
check(starts == 666666 and below == 666665 and cycled == [4] and capped == 0, "control S4 reproduction")
# the m=4 start: certify by hand through the E_- path 4 -> 2 -> 1 reversed: 1 -> 2 -> 5? no: 4 is reached from 8 (halving) or from ... 3n-1=4 impossible
print("  m = 4: G_- cycles on {4, 11}; but 1 -> 2 -> 4 is not an E_- path (2 -> 4 is not an arrow); 4 is reached from 8,")
print("  and 8 <- 16 <- 32 <- 11 <- 4 ... the E_- cycle 4->11->32->16->8->4; 1 reaches 4 via 1->2->5->14->7->20->10->29->86->43->128->64->32->16->8->4")
p4 = [1, 2, 5, 14, 7, 20, 10, 29, 86, 43, 128, 64, 32, 16, 8, 4]
for i in range(len(p4) - 1):
    check(is_E_minus_arrow(p4[i], p4[i + 1]), "path to 4")
print("  (verified, %d arrows). Hence Q2_- holds for all m <= %d, 3 !| m: FINITE-EXACT." % (len(p4) - 1, X3))
print("VERDICT: Q1_- and Q2_- both hold to 10^6, so the one-giant-SCC picture of E is SHEET-BLIND: E_- has it too,")
print("  while T_- has three cycles; the cycles are escaped through even->3n-1 arrows (the path lists above).")
tlog("S3 done")

# ---------------------------------------------------------------------------
# S4  (iv) carry statistics
# ---------------------------------------------------------------------------
hdr("S4  (iv) carry B_L: 2^K T^L(n) = 3^L n + B_L; B_L mod 4, mod 3^L, v_2(B_L)")
print("PROVED: B_L = b * (3^(L-1) + 3^(L-2) 2^(K_1) + ... + 2^(K_(L-1))) depends only on (b, word);")
print("  B_L is odd (first term odd, others even), so v_2(B_L) = 0 always, on both sheets;")
print("  B_L(word, -1) = -B_L(word, +1) exactly; every word is realised on both sheets (residue bijection, control S2).")
print()
print("cycles (word repeated m times, m = 1..3):")
print("sheet cycle-min | m | L K | B_L | B_L mod 4 | B_L mod 3^L | v_2(B_L) | gate B_L/(2^K-3^L)")
for b, cyc in [(1, cycle_plus[0])] + [(-1, c) for c in cycles_minus]:
    w1 = word_of(cyc[0], b, len(cyc))
    for m in range(1, 4):
        w = w1 * m
        B, K = carry(w, b)
        L = len(w)
        gate = Fraction(B, 2 ** K - 3 ** L)
        check(gate == cyc[0], "gate = cycle minimum")
        print("%+d %3d | %d | %2d %2d | %d | %d | %d | %d | %s" % (b, cyc[0], m, L, K, B, B % 4, B % 3 ** L, v2(B), gate))
# same word on the two sheets, random-looking starts
print()
print("same word, both sheets (word = plus-sheet word of n, L = 8):")
print("n | word | B_8(+) | B_8(-) | sum")
for n in [3, 7, 27, 97, 871, 6171, 77031, 837799]:
    w = word_of(n, 1, 8)
    Bp, K = carry(w, 1)
    Bm, _ = carry(w, -1)
    check(Bp + Bm == 0, "sign mirror")
    print("%d | %s | %d | %d | %d" % (n, w, Bp, Bm, Bp + Bm))
# distribution of B_L mod 4 and mod 3 over each sheet's OWN orbit words, n odd <= 2*10^5, L = 10
print()
print("distribution over odd n <= 200001 of the sheet's own orbit word (L = 10):")
LL = 10
for b in (1, -1):
    d4 = {}
    d3 = {}
    for n in range(1, 200002, 2):
        w = word_of(n, b, LL)
        B, K = carry(w, b)
        d4[B % 4] = d4.get(B % 4, 0) + 1
        d3[B % 3] = d3.get(B % 3, 0) + 1
    print("  sheet %+d: B mod 4 -> %s ; B mod 3 -> %s" % (b, sorted(d4.items()), sorted(d3.items())))
print("  (mod 4 the two sheets are swapped by 1 <-> 3, mod 3 by 1 <-> 2: the mirror B -> -B; the counts are the")
print("   counts of words with the property, identical on both sheets by the residue bijection up to the sign flip)")
# strict check of the mirror at the level of counts
c4 = {}
for b in (1, -1):
    d = {}
    for n in range(1, 200002, 2):
        w = word_of(n, b, LL)
        B, _ = carry(w, b)
        d[(b * B) % 4] = d.get((b * B) % 4, 0) + 1
    c4[b] = sorted(d.items())
print("  after multiplying by b (i.e. |B| mod 4): plus %s, minus %s" % (c4[1], c4[-1]))
print("  these differ only because the orbits of n <= 200001 realise the words with slightly different multiplicities")
print("  (a word of length 10 with exponent sum K is a class of odd n mod 2^(K+1); [1,200001] cuts these classes unevenly).")
print("  exact check on COMPLETE classes: odd n < 2^19 whose word of length 10 has K_10 <= 18 (audit A4 recomputes it):")
cc = {}
for b in (1, -1):
    c4 = {}
    c3 = {}
    for n in range(1, 1 << 19, 2):
        w = word_of(n, b, LL)
        B, K = carry(w, b)
        if K <= 18:
            a_ = abs(B)
            c4[a_ % 4] = c4.get(a_ % 4, 0) + 1
            c3[a_ % 3] = c3.get(a_ % 3, 0) + 1
    cc[b] = (sorted(c4.items()), sorted(c3.items()))
    print("    sheet %+d: |B| mod 4 -> %s ; |B| mod 3 -> %s" % (b, cc[b][0], cc[b][1]))
check(cc[1] == cc[-1], "exact mirror of |B| distributions on complete residue classes")
print("    identical on both sheets, exactly, as the word bijection requires.")
print("VERDICT (iv): B_L mod 4, mod 3^L and v_2(B_L) are word functions times b; sign-specific content = sign(B_L) = b only.")
tlog("S4 done")

# ---------------------------------------------------------------------------
# S5  (v) greedy words versus cycle words
# ---------------------------------------------------------------------------
hdr("S5  (v) greedy inverse G_b(m) = (2^k m - b)/3 versus the cycle words")
print("greedy k by residue m mod 9 (k minimal >= 0 with integer result not 0 mod 3):")
for b in (1, -1):
    tab = []
    for r in range(9):
        if r % 3 == 0:
            continue
        tab.append((r, greedy_k(r + 9, b)))
    print("  b=%+d: %s" % (b, tab))
check([greedy_k(r + 9, -1) for r in (1, 2, 4, 5, 7, 8)] == [1, 0, 3, 0, 1, 2], "control's k_- table")
print()
print("cycle words (halving exponents along the odd skeleton) and the greedy reverse:")
for b, cyc in [(1, cycle_plus[0])] + [(-1, c) for c in cycles_minus]:
    L = len(cyc)
    w = word_of(cyc[0], b, L)
    rev_word = tuple(reversed(w))
    # greedy reversibility: at each node cyc[i+1], is the greedy predecessor cyc[i] with k = w[i]?
    rows = []
    allg = True
    for i in range(L):
        tgt = cyc[(i + 1) % L]
        pred, k = G(tgt, b)
        ok = (pred == cyc[i] and k == w[i])
        allg = allg and ok
        rows.append("%d<-%d(k=%d)%s" % (tgt, pred, k, "" if ok else "!=cycle pred %d(k=%d)" % (cyc[i], w[i])))
    print("  sheet %+d cycle %s: word %s, K/L = %d/%d, reverse word %s, greedy-reversible: %s" % (b, cyc, w, sum(w), L, rev_word, allg))
    print("    greedy predecessors: " + "; ".join(rows))
print()
print("greedy orbits (12 steps) from the cycle minima, value(k):")
for b, start in [(1, 1), (-1, 1), (-1, 5), (-1, 17), (-1, 7), (-1, 91)]:
    x = start
    seq = []
    for _ in range(12):
        x, k = G(x, b)
        seq.append("%d(%d)" % (x, k))
    print("  sheet %+d from %d: %s" % (b, start, " ".join(seq)))
# 2-cycles of the full inverse relation m -> (2^k m - b)/3 (any k >= 0): m (2^(a+c) - 9) = b (2^c + 3)
print()
print("2-cycles of the inverse relation (any k >= 0): m(2^(a+c) - 9) = b(2^c + 3), a+c <= 40:")
for b in (1, -1):
    sols = []
    for s in range(0, 41):
        for c in range(0, s + 1):
            a = s - c
            den = 2 ** s - 9
            num = b * (2 ** c + 3)
            if den != 0 and num % den == 0 and num // den > 0:
                m = num // den
                m2 = (2 ** a * m - b)
                if m2 % 3 == 0:
                    m2 //= 3
                    if m2 > 0 and (2 ** c * m2 - b) == 3 * m:
                        sols.append((m, m2, a, c))
    print("  b=%+d: %s" % (b, sols))
    if b == 1:
        check(sols == [(1, 1, 2, 2)], "plus 2-cycles")
    else:
        check(sorted(set((min(x, y), max(x, y)) for x, y, _, _ in sols)) == [(1, 1), (4, 11), (5, 7)], "minus 2-cycles")
print("  the sign gate: plus needs 2^(a+c) > 9, minus 2^(a+c) < 9 (a+c <= 3); the minus 2-cycles are {1},{4,11},{5,7}")
print("  plus census complete for ALL a, c (PROVED): m >= 1 needs 2^c (2^a - 1) <= 12, so c <= 3 when a >= 1 (then 2^(a+c) <= 20,")
print("  covered); a = 0 needs 2^c - 9 | 2^c + 3, i.e. 2^c - 9 | 12, and 2^c in {10,11,12,13,15,21} is impossible.")
print("  and {4,11} (c = 0 i.e. an even->3n-1 arrow) is the G_- 2-cycle of S3; {5,7} is the T_- cycle, not greedy.")
# word statistics of greedy words: mean k over the first J greedy steps, both sheets, starts <= 10^5
print()
print("greedy word statistics: mean k of the first step over 3 !| m <= 10^5, and mean over 6 steps:")
for b in (1, -1):
    s1 = 0
    s6 = 0
    cnt_ = 0
    for m in range(1, 100001):
        if m % 3 == 0:
            continue
        cnt_ += 1
        x = m
        for j in range(6):
            x, k = G(x, b)
            if j == 0:
                s1 += k
            s6 += k
    print("  b=%+d: starts %d, mean k_1 = %.6f, mean k over 6 steps = %.6f" % (b, cnt_, s1 / cnt_, s6 / (6 * cnt_)))
rc9 = {}
for m in range(1, 100001):
    if m % 3:
        rc9[m % 9] = rc9.get(m % 9, 0) + 1
kp9 = {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}
km9 = {1: 1, 2: 0, 4: 3, 5: 0, 7: 1, 8: 2}
mp9 = Fraction(sum(kp9[r] * rc9[r] for r in kp9), sum(rc9.values()))
mm9 = Fraction(sum(km9[r] * rc9[r] for r in km9), sum(rc9.values()))
print("  the first-step means are exact rationals from the residue counts mod 9 of [1,10^5] (%s):" % sorted(rc9.items()))
print("    plus %s, minus %s, difference %s: the residue 1 occurs once more than the others and k_+(1) = 2, k_-(1) = 1"
      % (mp9, mm9, mp9 - mm9))
check(mp9 == Fraction(77779, 66667) and mm9 == Fraction(77778, 66667), "exact greedy means")
print("  (the six-step means differ in the fifth decimal for the same finite-range reason; the greedy residue chain")
print("   is conjugate under r -> -r mod 9, cited SCC Thm 6.2 / 7.1, so the full-period statistics coincide exactly)")
print("VERDICT (v): the greedy mechanism is sheet-blind; the minus cycles {5,7},{17..91} are not greedy-reversible,")
print("  the two trivial cycles are; the only sign-specific fact is the gate 2^(a+c) vs 9, i.e. the sign law.")
tlog("S5 done")

# ---------------------------------------------------------------------------
# S6  (vi) the word-function theorem and the residue-of-minimum law
# ---------------------------------------------------------------------------
hdr("S6  (vi) word-function invariants are blind; the order-only residue-of-minimum law")
print("PROVED (word-function theorem): if an invariant I(n) is a function of (b, parity word of n), then the")
print("  number of residues mod 2^J with I in any set is the same on both sheets after b -> -b (control S2 bijection).")
print("  Hence a sign-specific invariant must use the ORDER of Z (positivity or size), not residues or words alone.")
print()
print("residue-of-minimum law (order + residue): for a T_b cycle with minimum n_0 > 1 the step out of n_0 ascends,")
print("  (3n_0+b)/2^k > n_0 forces k = 1, so 3n_0+b = 2 mod 4: plus n_0 = 3 mod 4, minus n_0 = 1 mod 4;")
print("  the maximum M > 1 descends, forcing k >= 2, i.e. 3M+b = 0 mod 4: plus M = 1 mod 4, minus M = 3 mod 4.")
for b, cyc in [(-1, c) for c in cycles_minus[1:]]:
    n0 = min(cyc)
    M = max(cyc)
    _, k0 = T(n0, b)
    _, kM = T(M, b)
    print("  sheet %+d cycle %s: min %d = %d mod 4 (k=%d), max %d = %d mod 4 (k=%d)" % (b, cyc, n0, n0 % 4, k0, M, M % 4, kM))
    check(n0 % 4 == 1 and k0 == 1 and M % 4 == 3 and kM >= 2, "residue-of-minimum law")
print("  the plus fixed point 1 has k = 2 (3+1 = 4), the minus fixed point 1 has k = 1 (3-1 = 2): the root asymmetry again.")
print("  Decisive test (plus sheet, order-only): the last odd value > 1 before a T_+ orbit hits 1 is a descending")
print("  (k >= 2) node, hence 1 mod 4 by the law; indeed it is (4^j - 1)/3 = 1 + 4 + ... + 4^(j-1) = 1 mod 4 for the")
print("  j with 3n+1 = 4^j.  The census over odd 3 <= n <= 10^6 confirms it and shows the law has no exclusion content:")
lastodd = {}
for n in range(3, 10 ** 6 + 1, 2):
    x = n
    prev = x
    steps = 0
    while x != 1:
        prev = x
        x, _ = T(x, 1)
        steps += 1
        if steps > 10000:
            raise RuntimeError("plus orbit cap")
    lastodd[prev % 4] = lastodd.get(prev % 4, 0) + 1
print("  last odd value before 1, residue mod 4 histogram over odd 3 <= n <= 10^6:", sorted(lastodd.items()))
check(set(lastodd) == {1}, "last odd before 1 is 1 mod 4")
print("  every plus orbit to 10^6 ends through a 1 mod 4 node (exactly the maximum-type residue), as the law says;")
print("  a hypothetical plus cycle would also obey it with its own min = 3 mod 4 and max = 1 mod 4.")
print()
print("conjugation check: the law for the minus sheet is the law for negative n on the plus sheet (T_+(-n) = -T_-(n)),")
print("  and -1 mod 4 = 3 mod 4: the residue-of-minimum law is the sign law read mod 4.  It is SIGN-SPECIFIC in its")
print("  statement and reduces to the sign law; it is not new.")
print()
print("SURVIVING SIGN-SPECIFIC CONTENT (this lane's conclusion): positivity of n on the chosen sheet, equivalently")
print("  sign(B_L) = b, equivalently which side of log_2 3 a cycle's K/L must lie on; plus the finite verification")
print("  floor (2^68 plus, 10^7 minus, cited).  Every candidate (i)-(vi) is either a word function (blind), a")
print("  conjugation-image of the sign law, or a finite fact.  The ONE root asymmetry is S1: at n = 1 and only there")
print("  the two sheets' odd arrows fail strictness differently (halved: diagonal (1,1) vs companion 0;")
print("  unhalved: strict 1+3 vs diagonal 1+1), and 1 is fixed by the plus companion map (n+1)/2, killed by (n-1)/2.")
tlog("S6 done")

print()
print("ALL CHECKS PASSED")
