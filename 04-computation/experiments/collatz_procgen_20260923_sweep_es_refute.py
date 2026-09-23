#!/usr/bin/env python3
"""collatz_procgen_20260923_sweep_es_refute.py  (HYP-9121)

E_S with S = {6 mod 8}: x -> 3x+1 (x odd, forced); x -> x/2 (x even); and, only at even x = 6 mod 8,
the excursion x -> 3x+1 -> 9x+4.

(1) Lemma F: the forward E_S-orbit F of -1 is a finite set of 16 negative integers (closed under E_S).
(2) Lemma G: every cycle of the E_S graph on F has weight 3(#M) - 2(#H) >= 0; hence every E_S path from
    -1 has 3a >= 2b - C_F (C_F from shortest-path potentials); alpha_S(b) = least #M at the b-th halving.
(3) Theorem P_S family: x_c = -1 - 2^i c/3^j with j <= alpha_S(i-3) and 2^i c/3^j < (1-3^-j)/2 is
    E_S-hostile; exponentially many c at each scale.  Lower bound for |Bad_m(E_S)|.
(4) Cross-checks: membership of the family in the independently computed Bad_55 (thread DFS dump),
    and an exact E_S hostility prover (Lemma S / Lemma C of the dimension lane) on sample points.
usage: python3 ..._es_refute.py [BAD55_DUMP]
"""
import sys, math, collections
from fractions import Fraction as Fr

def es_succ(v):
    """E_S successor edges from an integer v: list of (w, 'M'/'H')"""
    out = []
    if v % 2 == 0: out.append((v // 2, 'H'))
    if v % 2 != 0 or v % 8 == 6: out.append((3 * v + 1, 'M'))
    return out

# (1) closure
F = {-1}; stack = [-1]
while stack:
    v = stack.pop()
    for w, _ in es_succ(v):
        if w not in F:
            F.add(w); stack.append(w)
            if len(F) > 10**6: raise SystemExit("not finite")
F = sorted(F, reverse=True)
print("Lemma F: forward E_S-orbit of -1 =", F, "(%d values, closed)" % len(F))
assert all(w in F for v in F for w, _ in es_succ(v))

# (2) weights: M -> +3, H -> -2 ; Bellman-Ford for negative cycles and potentials
idx = {v: k for k, v in enumerate(F)}
edges = [(idx[v], idx[w], 3 if t == 'M' else -2) for v in F for w, t in es_succ(v)]
n = len(F); INF = float('inf')
# all-pairs shortest paths (Floyd); negative cycle iff dist[u][u] < 0
d = [[INF] * n for _ in range(n)]
for u in range(n): d[u][u] = 0
for u, w, c in edges: d[u][w] = min(d[u][w], c)
for k in range(n):
    for u in range(n):
        if d[u][k] == INF: continue
        for w in range(n):
            if d[k][w] == INF: continue
            if d[u][k] + d[k][w] < d[u][w]: d[u][w] = d[u][k] + d[k][w]
neg = [F[u] for u in range(n) if d[u][u] < 0]
print("Lemma G: negative-weight cycles:", neg if neg else "none")
CF = -min(d[idx[-1]][w] for w in range(n) if d[idx[-1]][w] < INF)
print("  C_F = max over E_S paths from -1 of (2b - 3a) =", CF, " => a >= (2b - C_F)/3")
# simple cycles by DFS (small graph)
cyc = set()
def dfs(start, v, path, wsum, am, ah):
    for w, t in es_succ(v):
        if w == start:
            cyc.add((tuple(sorted(path)), am + (t == 'M'), ah + (t == 'H'))); continue
        if w in path or idx[w] < idx[start]: continue
        dfs(start, w, path + [w], 0, am + (t == 'M'), ah + (t == 'H'))
for s in F: dfs(s, s, [s], 0, 0, 0)
print("  simple cycles (vertex set size, #M, #H, rate #M/#H):",
      sorted(set((len(p), a, h, round(a / h, 4)) for p, a, h in cyc)))

# exact alpha_S(b) by DP on F (states after a halving; the start -1 counts as b = 0)
def alpha_S(bmax):
    cur = {-1: 0}; out = []
    for b in range(1, bmax + 1):
        nxt = {}
        for v, a in cur.items():
            # from v: sequences of M-edges (only M at odd or 6 mod 8) then one H
            frontier = [(v, a)]; seen = set()
            while frontier:
                u, au = frontier.pop()
                for w, t in es_succ(u):
                    if t == 'H':
                        if w not in nxt or nxt[w] > au: nxt[w] = au
                    else:
                        if (w, au + 1) not in seen:
                            seen.add((w, au + 1)); frontier.append((w, au + 1))
        cur = nxt; out.append(min(cur.values()))
    return out
AS = alpha_S(400)
th = math.log(2) / math.log(3)
print("alpha_S(b) - ceil((b+1) log3 2), b=1..80:", [AS[b-1] - math.ceil((b+1) * th) for b in range(1, 81)])
print("check a >= (2b - C_F)/3 for b <= 400:", all(3 * AS[b-1] >= 2 * b - CF for b in range(1, 401)))
print("alpha_S(b) - 2b/3 at b = 100, 200, 300, 400:", [round(AS[b-1] - 2 * b / 3, 3) for b in (100, 200, 300, 400)])


# (2b) loops through -1 in E_S: every closed walk through -1 uses -1 ->M -2 and -2 ->H -1 (weight +1), so
#      3a >= 2K + 1.  Exact least a for each K by DP on the finite graph F (all such loops live in F).
def min_loop_a(Kmax):
    INFa = 10**9
    # state: (vertex, halvings so far) -> least multiplications; start at -1 with 0,0; count halvings
    best = {(-1, 0): 0}; res = {}
    order = sorted(F)
    for K in range(0, Kmax + 1):
        # relax M-edges (no halving) to a fixed point at this K level (graph of M-edges is acyclic on F? use iteration)
        changed = True
        while changed:
            changed = False
            for v in F:
                a = best.get((v, K))
                if a is None: continue
                for w, t in es_succ(v):
                    if t == 'M' and best.get((w, K), INFa) > a + 1:
                        best[(w, K)] = a + 1; changed = True
        for v in F:
            a = best.get((v, K))
            if a is None: continue
            for w, t in es_succ(v):
                if t == 'H' and best.get((w, K + 1), INFa) > a:
                    best[(w, K + 1)] = a
        if K >= 1 and (-1, K) in best: res[K] = best[(-1, K)]
    return res
ml = min_loop_a(80)
a0 = lambda K: math.ceil((K + 1) * th)
exist = [K for K in range(1, 81) if ml.get(K) == a0(K)]
print("E_S loops through -1: K with a loop of a0(K) multiplications (K<=80):", exist)
print("  check 3a >= 2K+1 on every loop (K<=80):", all(3 * ml[K] >= 2 * K + 1 for K in ml))
print("  least a minus a0(K), K=1..80:", [ml.get(K, None) - a0(K) if K in ml else None for K in range(1, 81)])

# (3) the family and the count bound
def family(i):
    j = AS[i - 4]            # alpha_S(i-3)
    cmax = 0
    # largest c with 2^i c / 3^j < (1 - 3^-j)/2  <=>  2^(i+1) c 3^j < 3^j (3^j - 1)  <=> 2^(i+1) c < 3^j - 1
    cmax = (3**j - 2) // 2**(i + 1)
    return j, cmax
print("scale i: j = alpha_S(i-3), number of hostile perturbations c_max(i):")
for i in (10, 20, 40, 60, 80, 100, 150, 200):
    j, cm = family(i); print("  i=%d j=%d c_max=%d  (log2 c_max = %.2f)" % (i, j, cm, math.log2(cm) if cm > 0 else float('-inf')))
def lower_bound(m):
    best = 1
    for i in range(4, m):
        j, cm = family(i)
        best = max(best, min(cm, 2**(m - i)))
    return best
print("|Bad_m(E_S)| >= :", [(m, lower_bound(m)) for m in (40, 55, 64, 100, 200, 300, 400)])
print("  exponent log2(bound)/m at m = 200, 300, 400:", [round(math.log2(lower_bound(m)) / m, 4) for m in (200, 300, 400)])

# (3b) explicit bound check: exact lower bound vs 2^(0.0536 m - 6), 200 <= m <= 2000
AS2 = alpha_S(2100)
print("  alpha_S(b) >= ceil(2b/3) for b <= 2100:", all(AS2[b-1] >= -(-2*b//3) for b in range(1, 2101)))
def cmax2(i):
    j = AS2[i - 4]; return (3**j - 2) // 2**(i + 1)
worst = None
for m in range(200, 2001):
    best = 1
    for i in range(4, m):
        v = min(cmax2(i), 2**(m - i))
        if v > best: best = v
    r = math.log2(best) - (0.0536 * m - 6)
    if worst is None or r < worst[0]: worst = (r, m)
print("  min_{200<=m<=2000} [log2(bound) - (0.0536 m - 6)] = %.3f at m = %d" % worst)

# (4a) cross-check with the thread-DFS dump Bad_55 (independent code path)
if len(sys.argv) > 1:
    M = 55; bad = set(int(l) for l in open(sys.argv[1]) if l.strip())
    tot = hit = 0
    for i in range(4, 52):
        j, cm = family(i)
        for c in range(1, min(cm, 2**(M - i)) + 1):
            x = Fr(-1) - Fr(2**i * c, 3**j)
            cls = (x.numerator * pow(x.denominator, -1, 2**M)) % 2**M
            tot += 1; hit += cls in bad
    print("(4a) family members (i<52, c<=min(c_max,2^(55-i))): %d; in Bad_55 of the thread DFS: %d" % (tot, hit))

# (4b) exact E_S hostility prover (dimension-lane Lemmas S and C), sample points
def es_hostile(x, cap=2 * 10**6):
    """DFS over E_S paths from the negative rational x (odd denominator); SAFE at an integer v <= -1 with
    R > |v| (or v odd and R > |v| - 1/3); DESCENT if 3^a < 2^b at a halving; ESCAPE if a value >= 0."""
    stack = [(x, 0, 0)]; nodes = 0
    while stack:
        v, a, b = stack.pop(); nodes += 1
        if nodes > cap: return 'undecided'
        if v >= 0: return 'escape'
        if v.denominator == 1:
            R = Fr(3**a, 2**b)
            if R > -v or (v.numerator % 2 and R > -v - Fr(1, 3)): continue
        num = v.numerator
        if num % 2 == 0:
            w = v / 2
            if 3**a < 2**(b + 1): return 'descent'
            stack.append((w, a, b + 1))
            if (num * pow(v.denominator, -1, 8)) % 8 == 6:
                stack.append((3 * v + 1, a + 1, b))
        else:
            stack.append((3 * v + 1, a + 1, b))
    return 'hostile'
samples = []
for i in (10, 60):                     # the scales i < 80 with c_max(i) >= 1 whose path trees fit the prover's cap
    j, cm = family(i)
    for c in range(1, cm + 1): samples.append((i, j, c))
res = collections.Counter()
for i, j, c in samples:
    x = Fr(-1) - Fr(2**i * c, 3**j); r = es_hostile(x, cap=10**7); res[r] += 1
    if r != 'hostile': print("  sample", i, j, c, r)
print("(4b) exact prover on %d sample family members:" % len(samples), dict(res))

# (5) cross-check of (2b) with the tight DP of ..._sweep_tightdp.c (ES_FILTER=1), if the binary is given
import os, subprocess
exe = os.environ.get('TIGHTDP')
if exe and os.path.exists(exe):
    from fractions import Fraction
    agree = 0; tested = 0
    for K in range(10, 81):
        n = K + 1; a = 0
        while 3**a <= 2**n: a += 1
        G = Fraction(3**a + 1, 2**n) - 1
        if G >= Fraction(1, 2): continue          # tightness lemma T3(a) needs eta < log_3(3/2)
        env = dict(os.environ, ES_FILTER='1')
        r = subprocess.run([exe, '-1', '1', str(a), str(K), '%.17e' % float(G), str(10**7), '0', '0', '.'],
                           capture_output=True, text=True, env=env)
        found = 'found=1' in r.stdout; tested += 1
        agree += (found == (K in exist))
    print("(5) tight DP with ES_FILTER=1 vs exact F-DP on loops through -1: %d/%d agree" % (agree, tested))

# (6) other partial-choice sets: is the forward E_S-orbit of -1 finite (cap 10^4 vs 10^6)?
def orbit_size(Sres, J, CAP):
    seen = {-1}; dq = collections.deque([-1])
    while dq:
        v = dq.popleft(); nb = []
        if v % 2 == 0: nb.append(v // 2)
        if v % 2 != 0 or (v % (1 << J)) in Sres: nb.append(3 * v + 1)
        for w in nb:
            if abs(w) <= CAP and w not in seen: seen.add(w); dq.append(w)
    return len(seen)
print("(6) size of the forward orbit of -1 (values <= 10^4 / <= 10^6):")
for name, S, J in [("6 mod 8", {6}, 3), ("6,14 mod 16", {6, 14}, 4), ("2 mod 8", {2}, 3), ("0 mod 8", {0}, 3),
                   ("0 mod 4", {0}, 2), ("54 mod 64", {54}, 6), ("2 mod 4", {2}, 2), ("all evens (E)", {0}, 1)]:
    print("   S = %-14s %7d %7d" % (name, orbit_size(S, J, 10**4), orbit_size(S, J, 10**6)))
