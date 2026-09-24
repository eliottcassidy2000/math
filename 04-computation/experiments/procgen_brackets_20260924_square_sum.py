#!/usr/bin/env python3
"""Part 3: the square-sum problem, its small-n -> large-n transition, and a mechanism test against
Collatz-type small-number phenomena.

G_n: vertices 1..n, edge {a,b} (a != b) iff a+b is a square.  A square-sum chain = Hamiltonian path of G_n;
a square-sum loop = Hamiltonian cycle.  Target variants: triangular numbers, squares+doubled squares, odd squares.

Literature (read this session: OEIS A090461, A071983, A071984, A090460 entries; the reversed-square-sum README):
  paths exist exactly for n in {1, 15, 16, 17, 23} and every n >= 25 (claimed PROVED by R. Gerbicz, Mersenneforum,
  Jan 2018, elementary induction n -> 49n + r plus finite tables; forum thread login-gated, not read directly);
  loops exist exactly for n >= 32 (same proof; n <= 30 excluded by a vertex of degree <= 1, n = 31 by hand).

Sections
 (a) exact existence for n <= 60 and obstruction types            FINITE-EXACT
 (b) Hamiltonian path / cycle counts n = 15..38 against A071983/A071984  FINITE-EXACT (independent code)
 (c) transition metrics: degrees, leaves, count growth               FINITE-EXACT
 (d) the self-similar engine: 49*a +- c chains partition [25, 49n+24]  PROVED identity, checked
 (e) mechanism test 1: other target sets (triangular numbers, squares+doubled squares, cubes from the literature):
     the loop threshold sits where the typical degree c*n^alpha reaches ~2.3
 (f) mechanism test 2: cycles of (qn +- 1)/2 are governed by |2^K - q^a|, not by any smooth density
Runtime ~ 1-3 minutes.  Memory small.
"""
import os, sys, subprocess, tempfile, math, time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="procgen_brackets_ss_")
SS = os.path.join(TMP, "squaresum"); PB = os.path.join(TMP, "pairings")
subprocess.run(["cc", "-O2", "-o", SS, os.path.join(HERE, "procgen_brackets_20260924_squaresum.c")], check=True)
subprocess.run(["cc", "-O2", "-o", PB, os.path.join(HERE, "procgen_brackets_20260924_pairings.c"), "-lm"], check=True)
def run(binp, *args, timeout=None):
    return subprocess.run([binp] + [str(a) for a in args], capture_output=True, text=True, check=True, timeout=timeout).stdout.splitlines()

# OEIS data (fetched 2026-09-24): A071983 (undirected Hamiltonian paths, n = 15..), A071984 (loops, n = 32..)
A071983 = [1,1,1,0,0,0,0,0,3,0,10,12,35,52,19,20,349,392,669,4041,17175,12960,14026,11889,29123,39550]
A071984 = [1,1,11,57,31,20,25,50,64,464]
OEIS_PATH_N = {1, 15, 16, 17, 23} | set(range(25, 200))

print("=" * 100)
print("PART 3. The square-sum problem and the microcosm/macrocosm comparison")
print("=" * 100)

# ------------------------------------------------------------------------------------------ (a)
print("\n(a) Exact existence (exhaustive DFS with dead-end pruning; obstruction reported when trivial)")
t0 = time.time()
rows = {}
for l in run(SS, "exist", 1, 60, 0, 10**9):
    f = l.split()
    n = int(f[4]); rows[n] = dict(mindeg=int(f[6]), leaves=int(f[8]), iso=int(f[10]), path=f[12], cycle=f[14])
paths = sorted(n for n, r in rows.items() if r["path"] == "YES")
cycles = sorted(n for n, r in rows.items() if r["cycle"] == "YES")
unk = sorted(n for n, r in rows.items() if "UNKNOWN" in (r["path"], r["cycle"]))
print(f"   n <= 60 with a square-sum chain: {paths[:6]} + all of [25..60]: {set(range(25, 61)) <= set(paths)}; matches A090461: {set(paths) == OEIS_PATH_N & set(range(1, 61))}")
print(f"   n <= 60 with a square-sum loop : first {cycles[:3]}, all of [32..60]: {cycles == list(range(32, 61))}; undecided: {unk}")
print("   obstruction types for the failures:")
for n in range(2, 32):
    r = rows[n]
    if r["path"] != "YES" or r["cycle"] != "YES":
        print(f"     n={n:2d}  mindeg {r['mindeg']} leaves {r['leaves']} isolated {r['iso']}   path {r['path']:<13} loop {r['cycle']}")
print(f"   ({time.time()-t0:.1f}s)")

# ------------------------------------------------------------------------------------------ (b)
print("\n(b) Counts of chains and loops against OEIS (independent exhaustive enumeration)")
t0 = time.time()
cnt = {}
for l in run(SS, "count", 15, 38, 0, 5 * 10**10, timeout=3000):
    f = l.split(); cnt[int(f[4])] = (int(f[6]), int(f[8]), f[9])
okp = all(cnt[n][0] == A071983[n - 15] for n in cnt if n - 15 < len(A071983))
okc = all(cnt[n][1] == A071984[n - 32] for n in cnt if n >= 32 and n - 32 < len(A071984))
print(f"   n = 15..38 chains: {[cnt[n][0] for n in sorted(cnt)]}")
print(f"   n = 32..38 loops : {[cnt[n][1] for n in sorted(cnt) if n >= 32]}")
print(f"   agree with A071983: {okp};  agree with A071984: {okc};  all exact: {all(v[2] == 'exact' for v in cnt.values())}  ({time.time()-t0:.1f}s)")

# ------------------------------------------------------------------------------------------ (c)
print("\n(c) Transition metrics")
def degs(n, target=lambda s: math.isqrt(s) ** 2 == s):
    d = [0] * (n + 1)
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            if target(a + b): d[a] += 1; d[b] += 1
    return d[1:]
last_leaf = max(n for n in range(2, 400) if min(degs(n)) <= 1)
print(f"   last n with a vertex of degree <= 1: {last_leaf}  (vertex 18 has the single partner 7 until 31 joins at n = 31)")
for n in (10, 25, 50, 100, 200, 400):
    d = degs(n)
    print(f"     n={n:4d}: min degree {min(d)}, mean degree {sum(d)/n:.2f}, (sqrt2-1)*sqrt(n) = {(math.sqrt(2)-1)*math.sqrt(n):.2f}, sqrt(2n)-1 = {math.sqrt(2*n)-1:.2f}")
g = [(n, cnt[n][0]) for n in sorted(cnt) if cnt[n][0] > 0]
print("   log(#chains)/n:", [(n, round(math.log(c) / n, 3)) for n, c in g if n >= 25][::3])
# the exact count of squares in dyadic windows, the common mechanism with the brackets
sq_in = lambda x: sum(1 for k in range(1, 3 * math.isqrt(2 * x) + 3) if x < k * k <= 2 * x)
odd_sq_in = lambda x: sum(1 for k in range(1, 3 * math.isqrt(2 * x) + 3, 2) if x < k * k <= 2 * x)
no_odd = [x for x in range(1, 200) if odd_sq_in(x) == 0]
one_sq = [x for x in range(1, 200) if sq_in(x) <= 1]
print(f"   x < 200 with no odd square in (x, 2x]: {no_odd}  (= S_2 of part 1 plus x = 1)")
print(f"   x < 200 with at most one square in (x, 2x]: {one_sq}")

# ------------------------------------------------------------------------------------------ (d)
print("\n(d) The self-similar engine behind the proof (Gerbicz): if a + b = s^2 then (49a + c) + (49b - c) = (7s)^2")
def a_chain(n):
    # a square-sum chain of [1,n] by DFS (small n)
    adj = {a: [b for b in range(1, n + 1) if b != a and math.isqrt(a + b) ** 2 == a + b] for a in range(1, n + 1)}
    order = sorted(adj, key=lambda v: len(adj[v]))
    sys.setrecursionlimit(10000)
    def rec(path, used):
        if len(path) == n: return list(path)
        for w in sorted(adj[path[-1]], key=lambda v: len(adj[v])):
            if w not in used:
                used.add(w); path.append(w)
                r = rec(path, used)
                if r: return r
                path.pop(); used.discard(w)
        return None
    for s in order:
        r = rec([s], {s})
        if r: return r
P = a_chain(25)
chains = [[49 * a + (c if p % 2 == 0 else -c) for p, a in enumerate(P)] for c in range(-24, 25)]
allv = sorted(v for ch in chains for v in ch)
ok_sq = all(math.isqrt(x + y) ** 2 == x + y for ch in chains for x, y in zip(ch, ch[1:]))
print(f"   base chain of [1,25]: {P}")
print(f"   the 49 chains T(c) = (49 a_1 + c, 49 a_2 - c, ...), c = -24..24, are square-sum chains: {ok_sq};")
print(f"   together they partition [25, 49*25+24] = [25, 1249] exactly: {allv == list(range(25, 1250))}")
print("   The remaining work of the proof is to glue the 49 chains with the 24 small numbers (Gerbicz's 'nice pairs');")
print("   that step is not re-verified here.  This exact self-similarity is the macrocosm engine of the problem.")

# ------------------------------------------------------------------------------------------ (e)
print("\n(e) Mechanism test 1: other target sets.  A vertex near n has about c*n^alpha partners.  If the transition")
print("    is a degree threshold, the last exception n* should sit where c*n*^alpha reaches the same small value.")
targets = [(0, "squares", math.sqrt(2) - 1, 0.5), (2, "triangular k>=2", 2 - math.sqrt(2), 0.5),
           (3, "squares+2squares", 1 - 1 / math.sqrt(2) + math.sqrt(2) - 1, 0.5)]
for tgt, name, c, al in targets:
    t0 = time.time()
    out = run(SS, "exist", 2, 127, tgt, 2 * 10**7, timeout=3000)
    res = {}
    for l in out:
        f = l.split(); res[int(f[4])] = (f[12], f[14], int(f[8]), int(f[6]))
    no_p = [n for n in res if res[n][0].startswith("NO")]; unk_p = [n for n in res if res[n][0] == "UNKNOWN"]
    no_c = [n for n in res if res[n][1].startswith("NO")]; unk_c = [n for n in res if res[n][1] == "UNKNOWN"]
    lp, lc = max(no_p), max(no_c)
    deg_p = max(n for n in res if res[n][0] in ("NO(isolated)", "NO(leaves>2)"))
    deg_c = max(n for n in res if res[n][1] == "NO(deg<2)")
    print(f"   {name:17s} c = {c:.3f}: last n without chain {lp:3d} (last degree obstruction {deg_p:3d}); last n without loop"
          f" {lc:3d} (last degree obstruction {deg_c:3d});  c*sqrt(n) at the loop threshold n = {lc+1}: {c*math.sqrt(lc+1):.2f};"
          f" undecided chains {unk_p}, loops {unk_c}  ({time.time()-t0:.0f}s)")
print("   cubes (literature: Rivera's Puzzle 311 via OEIS A071984): first cubic loop at n = 473;")
print(f"     c = 2^(1/3) - 1 = {2**(1/3)-1:.3f}, alpha = 1/3: c*473^(1/3) = {(2**(1/3)-1)*473**(1/3):.2f}")
print("   odd squares: a + b = 1 (mod 8) splits G_n into the residue pairs {0,1},{2,7},{3,6},{4,5} mod 8, so no chain")
print("   exists for any n >= 8 (a congruence obstruction, not a density one): thinning the target by congruences")
print("   disconnects the graph, which is why the density test uses unions and triangular numbers instead.")
# ------------------------------------------------------------------------------------------ (f)
print("\n(f) Mechanism test 2: small cycles of the Collatz-type maps n -> (qn + s)/2 (n odd), n/2 (n even)")
print("    For a cycle with a odd steps and K steps: n |2^K - q^a| = |B| <= K 2^(K-1) (q/2)^a; big cycles need 2^K ~ q^a.")
def stats(q, s, m):
    v, a, K, mx = m, 0, 0, m
    while True:
        if v % 2: v = (q * v + s) // 2; a += 1
        else: v //= 2
        K += 1; mx = max(mx, v)
        if v == m: return a, K, mx
t0 = time.time()
table = []
for q in list(range(3, 40, 2)) + [181, 1093]:
    for s in (1, -1):
        l = run(PB, "qmap", q, s, 20000, 10**15)[0]
        tail = l.split("|")[1].split()
        cyc = [tuple(int(x) for x in t.split(":")) for t in tail[5:]]
        for mn, ln, mxv, basin in cyc:
            if mn == 0 or (mn == 1 and ln <= 3):
                continue
            a, K, mx = stats(q, s, mn)
            rel = abs(2**K - q**a) / max(2**K, q**a)
            table.append((q, s, mn, ln, a, K, mx, rel))
print("    nontrivial cycles with minimum <= 20000 (q odd <= 39, plus 181 and 1093):")
print("      q  sign  min   len  odd  |2^K - q^a|/max   max element   max/min")
for q, s, mn, ln, a, K, mx, rel in table:
    print(f"    {q:4d}  {'+' if s > 0 else '-'}  {mn:6d} {ln:4d} {a:4d}   {rel:12.3e}   {mx:12d}   {mx/mn:9.1f}")
print(f"    ({time.time()-t0:.1f}s)  Every nontrivial cycle sits at a near-coincidence 2^K ~ q^a; no smooth threshold in q exists")
print("    (q = 181: 2^15 - 181^2 = 7 gives cycles with elements up to 55296; q = 1093, a base-2 Wieferich prime, has none.)")
