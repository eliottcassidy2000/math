#!/usr/bin/env python3
"""Part 2: the Skolem/graceful view of the Collatz shortcut map and the 'pairing ladder'.

T(n) = n/2 (n even), (3n+1)/2 (n odd).  |T(n) - n| = ceil(n/2): the pair {2i-1, 2i} shares the length i.
Pairing maps F_eps (offset 0: pairs {2i-1,2i}; offset 1: pairs {2i,2i+1}): n goes UP by its length iff
(n odd) XOR eps_i, else DOWN.  eps = 0 gives T (offset 0) and the 3n-1 map U (offset 1).

Sections
 (A) pair-sum invariant, AM-GM drift, uniqueness of q = 3        (PROVED; exhaustive checks)
 (B) graceful form: two perfect difference systems; truncations  (PROVED + FINITE-EXACT)
 (C) symmetries of the family: shift+complement, Collatz <-> 3n-1 antipodal  (PROVED; checked)
 (D) periodic members mod 2^K: cycles, measure preservation, exceptional classes; the pair-0 obstruction
 (E) random members (Monte Carlo)
 (F) single flips of T and of U: fragility is a small-number phenomenon
 (G) the provable rung: landing => down; minimal flip density (CP-SAT window optimum)
 (H) a density-zero set of flips that makes an orbit diverge (PROVED construction)
Needs: cc, the C helper procgen_brackets_20260924_pairings.c (compiled into a temp dir); ortools for (G).
Session collatz-procgen-20260923/24, brackets lane.  Runtime ~ 3-6 minutes; memory < 400 MB.
"""
import os, sys, subprocess, tempfile, math, random, statistics, time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="procgen_brackets_")
BIN = os.path.join(TMP, "pairings")
subprocess.run(["cc", "-O2", "-o", BIN, os.path.join(HERE, "procgen_brackets_20260924_pairings.c"), "-lm"], check=True)
CAP = str(4 * 10**18)

def run(*args):
    return subprocess.run([BIN] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout.splitlines()

def parse_cycles(line):
    """'... | ncyc k nesc e cycles min:len:max:basin ...' -> (nesc, [(min,len,max,basin)])"""
    tail = line.split("|", 1)[1].split()
    nesc = int(tail[3])
    cyc = [tuple(int(x) for x in t.split(":")) for t in tail[5:]]
    return nesc, cyc

def attracting(cyc):
    return [c for c in cyc if c[3] > 0]

def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
def U(n): return n // 2 if n % 2 == 0 else (3 * n - 1) // 2
def Fq(q, n): return n // 2 if n % 2 == 0 else (q * n + 1) // 2

print("=" * 100)
print("PART 2. The Skolem/graceful view of T and the pairing ladder")
print("=" * 100)

# ------------------------------------------------------------------------------------------ (A)
print("\n(A) Pair-sum invariant and the AM-GM drift (PROVED; checks below are exhaustive in range)")
N = 10**6
okA = all(abs(T(n) - n) == (n + 1) // 2 for n in range(1, 2 * N + 1))
okS = all(T(2 * i - 1) + T(2 * i) == 4 * i - 1 == (2 * i - 1) + 2 * i for i in range(1, N + 1))
okU = all(U(2 * i) + U(2 * i + 1) == 4 * i + 1 and {U(2 * i), U(2 * i + 1)} == {i, 3 * i + 1} for i in range(1, N + 1))
from fractions import Fraction as Fr
okP = all(Fr((3 * i - 1) * i, (2 * i - 1) * 2 * i) == Fr(3, 4) + Fr(1, 8 * i - 4) for i in range(1, 5000))
selfT = [i for i in range(1, 10**5) if {T(2 * i - 1), T(2 * i)} == {2 * i - 1, 2 * i}]
selfU = [i for i in range(0, 10**5) if {U(2 * i), U(2 * i + 1)} == {2 * i, 2 * i + 1}]
sums = all(sum(T(n) for n in range(1, 2 * M + 1)) == M * (2 * M + 1) for M in (1, 2, 3, 10, 1000, 12345))
print(f"   |T(n)-n| = ceil(n/2) for n <= 2e6: {okA};  T(2i-1)+T(2i) = 4i-1 = (2i-1)+2i for i <= 1e6: {okS}")
print(f"   3n-1 map: U(2i)+U(2i+1) = 4i+1, image pair {{i, 3i+1}}, i <= 1e6: {okU}")
print(f"   product ratio T(2i-1)T(2i)/((2i-1)2i) = 3/4 + 1/(8i-4) (exact, i < 5000): {okP}")
print(f"   pairs mapped onto themselves: T: i in {selfT} (the cycle {{1,2}});  U: i in {selfU} (pair {{0,1}}, both fixed)")
print(f"   sum_(n<=2M) T(n) = sum_(n<=2M) n = M(2M+1): {sums}")
# uniqueness of q = 3 among (qn+1)/2 and (qn-1)/2
uniq = []
for q in range(1, 100, 2):
    for sgn in (+1, -1):
        # a sum-preserving pairing must match odd n with the even m = (q-2)n + sgn  (from ((q-2)n+sgn)/2 = m/2)
        X = 4000
        partners = [(q - 2) * n + sgn for n in range(1, X, 2)]
        evens_needed = set(range(2 if sgn > 0 else 0, X - 1, 2))
        ok = all(m >= 0 and m % 2 == 0 for m in partners) and set(p for p in partners if p < X - 1) >= evens_needed \
             and len(set(partners)) == len(partners)
        if ok:
            uniq.append((q, sgn))
print(f"   odd q <= 99, sheets +-1: sum-preserving perfect matchings (odd n <-> even (q-2)n +- 1) exist only for {uniq}")
print("   (two evens or two odds can never share a pair: for q >= 3 both displacements then have the same sign;")
print("    for q = 1 every displacement is <= 0; so a sum-preserving pair is odd n with the even m = (q-2)n +- 1)")
Dq = [(q, sum(Fq(q, n) - n for n in range(1, 2 * 500 + 1)), (q - 3) * 500**2 // 2) for q in (1, 3, 5, 7, 9)]
print(f"   displacement over [1,2M] for (qn+1)/2, M = 500: (q, sum, (q-3)M^2/2) = {Dq}")
print("   => for q != 3 the displacement grows like M^2, so no partition into zero-sum blocks of bounded diameter exists;")
print("      (for q >= 5 unbounded 'Hilbert hotel' partitions do exist, see below, so the diameter bound is necessary).")
print("   Every pairing map F_eps preserves every pair sum as well (one member goes up by i, the other down by i).")
# unbounded-diameter zero-sum partitions do exist for q = 5, 7 (greedy 'Hilbert hotel' construction)
for q in (5, 7):
    used = set(); X = 20000; stuck = 0; maxdiam = 0; blocks = 0
    disp = lambda n: ((q - 2) * n + 1) // 2 if n % 2 else -(n // 2)
    for n0 in range(1, X + 1):
        if n0 in used: continue
        if n0 % 2:
            s2 = 2 * disp(n0)
            if s2 not in used: blk = [n0, s2]
            else:
                e = n0 + 1; blk = None
                while e <= s2:
                    if e % 2 == 0 and e not in used and (s2 - e) > 0 and (s2 - e) not in used and s2 - e != e:
                        blk = [n0, e, s2 - e]; break
                    e += 1
                if blk is None: stuck += 1; blk = [n0]
        else:
            k = n0 + 1
            while True:
                if k % 2 and k not in used:
                    mp = 2 * disp(k) - n0
                    if mp > 0 and mp % 2 == 0 and mp not in used and mp not in (n0, k):
                        blk = [n0, k, mp]; break
                k += 2
        assert len(blk) == 1 or sum(disp(x) for x in blk) == 0
        used.update(blk); blocks += 1; maxdiam = max(maxdiam, max(blk) - min(blk))
    print(f"   q = {q}: a greedy partition into finite zero-sum blocks covers [1,{X}] ({blocks} blocks, stuck {stuck}),"
          f" max block diameter {maxdiam}: zero-sum blocks exist, but only with unbounded diameter")

# ------------------------------------------------------------------------------------------ (B)
print("\n(B) Graceful form: down-edges {i,2i} and up-edges {2i-1,3i-1} each realise every difference exactly once")
Mx = 10**6
down = Counter(2 * i - i for i in range(1, Mx + 1)); up = Counter((3 * i - 1) - (2 * i - 1) for i in range(1, Mx + 1))
print(f"   each difference 1..{Mx} exactly once in each class: {set(down.values()) == {1} and set(up.values()) == {1} and len(down) == len(up) == Mx}")
# Skolem-graceful: halving forest H_n on [1,2n]
hs = all(sorted(abs(2 * i - i) for i in range(1, n + 1)) == list(range(1, n + 1)) for n in (1, 5, 50, 777))
print(f"   halving forest H_n = {{ {{i,2i}} : i <= n }} on vertex set [1,2n], identity labels: p = 2n vertices labelled 1..2n,")
print(f"   q = n edges with differences exactly 1..n -> Skolem graceful in the sense of Lee-Shee (1991): {hs}")
print("   (a Skolem sequence of order n is the same object with H_n replaced by the matching nK_2)")
rosa = True
for n in range(1, 60):
    Nn = 2 * n + 1; Es = set()
    for t in range(Nn):
        for i in range(1, n + 1):
            a_, b_ = (i + t) % Nn, (2 * i + t) % Nn
            e_ = (min(a_, b_), max(a_, b_))
            if e_ in Es or a_ == b_: rosa = False
            Es.add(e_)
    rosa = rosa and len(Es) == Nn * (Nn - 1) // 2
print(f"   Rosa: the 2n+1 translates of H_n mod 2n+1 partition the edges of K_(2n+1), n < 60: {rosa}")
# truncation to [1,M]
for M in (100, 10**4, 10**6):
    mult = Counter()
    for n in range(2, M + 1):          # edge n -> T(n) inside [1,M]; the edge {1,2} is shared by 1->2 and 2->1
        if T(n) <= M:
            mult[(n + 1) // 2] += 1
    twice = sum(1 for d, c in mult.items() if c == 2); once = sum(1 for d, c in mult.items() if c == 1)
    print(f"   induced T-graph on [1,{M}]: differences used twice {twice} (all d <= {max(d for d,c in mult.items() if c == 2)}), once {once} (up to {max(mult)});"
          f" predicted twice for d <= (M+1)/3 = {(M + 1) // 3}")
# backward tree from 1 within [1, M]
for M in (100, 10**4, 10**6):
    inR = bytearray(M + 1); inR[1] = inR[2] = 1; stack = [2]; seen = 2
    while stack:
        v = stack.pop()
        for u in (2 * v, (2 * v - 1) // 3 if (2 * v - 1) % 3 == 0 and (2 * v - 1) // 3 % 2 == 1 else 0):
            if 1 <= u <= M and not inR[u] and T(u) == v:
                inR[u] = 1; seen += 1; stack.append(u)
    R = [n for n in range(1, M + 1) if inR[n]]
    mult = Counter((n + 1) // 2 for n in R if n != 1)           # the edge {1,2} counted once (from 2)
    q = len(R) - 1
    print(f"   backward tree R_M from 1 inside [1,{M}]: |R_M| = {len(R)} ({len(R)/M:.3f} of [1,M]), q = {q} edges, labels up to {max(R)};"
          f" differences twice {sum(1 for c in mult.values() if c == 2)}, once {sum(1 for c in mult.values() if c == 1)},"
          f" max difference {max(mult)} (graceful would need labels in [0,q] and differences exactly 1..q)")
# literal gracefulness: any interval [a,b] whose induced T-graph is a tree with differences exactly 1..b-a
lit = []
for a in range(0, 60):
    for b in range(a + 1, 200):
        V = range(a, b + 1)
        E = set()
        for n in V:
            if n >= 1 and a <= T(n) <= b and T(n) != n:
                E.add((min(n, T(n)), max(n, T(n))))
        if len(E) != b - a:
            continue
        diffs = sorted(y - x for x, y in E)
        if diffs == list(range(1, b - a + 1)):
            # connectivity
            par = {v: v for v in V}
            def fd(x):
                while par[x] != x:
                    par[x] = par[par[x]]; x = par[x]
                return x
            for x, y in E:
                par[fd(x)] = fd(y)
            if len({fd(v) for v in V}) == 1:
                lit.append((a, b))
print(f"   intervals [a,b], a < 60, b < 200, whose induced T-graph is a tree with differences exactly 1..b-a: {lit}")
print("   PROOF that only {1,2} occurs: difference 1 forces the edge {1,2} (the only edge of length 1), so a <= 1;")
print("   difference q = b-a forces an edge from 2q-1 or 2q, so b >= 2q-1, i.e. q <= a+1 <= 2.  The same argument")
print("   rules out every pairing map F_eps beyond 3 edges (length-1 edges sit in {0,1,2,3}).")

# ------------------------------------------------------------------------------------------ (C)
print("\n(C) Symmetries of the pairing family")
def Fgen(eps, off, n):
    i = (n + 1) // 2 if off == 0 else n // 2
    up = (n & 1) ^ eps(i)
    return n + i if up else n - i
rng = random.Random(7)
bits = [rng.randrange(2) for _ in range(10**5 + 10)]
e1 = lambda i: bits[i]; e1c = lambda i: 1 - bits[i]
okshift = all(Fgen(e1, 0, n) + 1 == Fgen(e1c, 1, n + 1) for n in range(0, 10**5))
okneg = all(T(-n) == -U(n) for n in range(1, 10**4))
print(f"   F^(offset 0)_eps(n) + 1 = F^(offset 1)_(1-eps)(n+1) for a random eps, n < 1e5: {okshift}")
print(f"   hence 3n-1 (offset 1, eps = 0) is conjugate by n -> n+1 to the ALL-SWAPPED Collatz pairing (offset 0, eps = 1):")
print(f"   Collatz and 3n-1 are antipodal corners of the cube {{0,1}}^N of choices.  Negation T(-n) = -U(n): {okneg}")
print("   (so the all-swapped Collatz pairing has the cycles {4,6,9} and {16,...,135}: the 3n-1 cycles shifted by one)")

# ------------------------------------------------------------------------------------------ (D)
print("\n(D) Periodic members eps_i = c(i mod 2^K): all 2^(2^K) choices, both offsets")
print("    PROVED obstruction: the length-0 pair {-1,0} (offset 0) / {0,1} (offset 1) makes one of the 2-adic points")
print("    -1, 0 (offset 0) or 0, 1 (offset 1) an UP-fixed point with multiplier 3/2 (which one is decided by c(0)).")
print("    Its 2-adic neighbours rise for as many steps as they share digits, so every periodic member has a nonempty")
print("    exceptional set at every level: no bounded-lookahead descent certificate exists for any periodic pairing.")
per_rows = {}
t0 = time.time()
for K, off, NN, M in [(0, 0, 10**6, 20), (0, 1, 10**6, 20), (1, 0, 10**6, 20), (1, 1, 10**6, 20), (2, 0, 10**6, 20), (2, 1, 10**6, 20),
                      (3, 0, 10**6, 20), (3, 1, 10**6, 20), (4, 0, 10**4, 12), (4, 1, 10**4, 12)]:
    rows = []
    for line in run("periodic", K, off, NN, CAP, M):
        f = line.split()
        mask, ones, mp, land, exc = int(f[6]), int(f[8]), int(f[10]), int(f[12]), int(f[14])
        nesc, cyc = parse_cycles(line)
        rows.append((mask, ones, mp, land, exc, nesc, cyc))
    per_rows[(K, off)] = rows
    ntree = sum(1 for r in rows if len(attracting(r[6])) == 1 and r[5] == 0)
    nmp = sum(r[2] for r in rows); nland = sum(r[3] for r in rows); nesc_maps = sum(1 for r in rows if r[5] > 0)
    ncyc_dist = Counter(len(attracting(r[6])) for r in rows)
    excs = sorted(r[4] for r in rows)
    col = [r for r in rows if r[0] == 0][0]
    bigmin = max((c[0] for r in rows for c in attracting(r[6])), default=0)
    print(f"   K={K} off={off}: {len(rows)} maps (N={NN}); measure-preserving {nmp}; landing=>down {nland}; trees to N {ntree}"
          f" ({100*ntree/len(rows):.1f}%); maps with escapes {nesc_maps}; #attracting cycles {dict(sorted(ncyc_dist.items()))};"
          f" largest cycle minimum {bigmin}; exceptional classes at 2^{M}: min {excs[0]} median {excs[len(excs)//2]} max {excs[-1]}"
          f" (eps=0 member: {col[4]})")
print(f"   ({time.time()-t0:.1f}s)")
# periodic members with escaping orbits: up-frequency under Haar measure (Monte Carlo on random 2-adic integers)
def upfreq(mask, K, off, trials=400, steps=300, seed=11):
    rr = random.Random(seed); ups = tot = 0
    for _ in range(trials):
        x = rr.getrandbits(steps + K + 40) | (1 << (steps + K + 39))
        for _ in range(steps):
            i = (x + 1) // 2 if off == 0 else x // 2
            up = (x & 1) ^ ((mask >> (i % (1 << K))) & 1)
            x = x + i if up else x - i
            ups += up; tot += 1
    return ups / tot
esc_maps = [r for r in per_rows[(3, 0)] if r[5] > 0]
print(f"   K=3 offset 0: the {len(esc_maps)} maps with orbits beyond 4e18 (divergence-like), with the up-frequency of")
print(f"   Haar-random 2-adic orbits (divergence needs > log2/log3 = {math.log(2)/math.log(3):.4f}):")
print("     " + ", ".join(f"mask {r[0]} (mp {r[2]}, {r[5]} n escape, up {upfreq(r[0], 3, 0):.3f})" for r in esc_maps[:12]))
fr = [upfreq(r[0], 3, 0) for r in per_rows[(3, 0)] if r[5] == 0][:40]
print(f"     maps without escapes (first 40): up-frequency range [{min(fr):.3f}, {max(fr):.3f}]")
# thinnest exceptional sets (K = 2, 3, offset 0): growth with the level
print("   exceptional-class growth for the thinnest periodic maps (offset 0, K = 3) and for Collatz:")
rows = per_rows[(3, 0)]
thin = sorted(rows, key=lambda r: r[4])[:4]
for r in [rows[0]] + thin:
    out = run("cert", r[0], 3, 0, 28)
    seq = [(int(l.split()[8]), l.split()[10]) for l in out]
    print(f"     mask {r[0]:3d} (ones {r[1]}, mp {r[2]}): level:count {seq}  cycles(attracting) {[(c[0], c[1]) for c in attracting(r[6])]}")

# ------------------------------------------------------------------------------------------ (E)
print("\n(E) Random members: eps_i i.i.d. Bernoulli(p) relative to Collatz (p = 1/2 is the uniform random pairing)")
t0 = time.time()
for p, seeds, NN in [(0.5, 300, 10**5), (0.5, 12, 10**6), (0.1, 300, 10**5), (0.01, 300, 10**5), (0.001, 300, 10**5)]:
    trees = 0; ncs = []; mins = []; maxs = []; lens = []; esc = 0
    for s in range(1, seeds + 1):
        line = run("random", s, p, 0, NN, CAP)[0]
        nesc, cyc = parse_cycles(line)
        at = attracting(cyc)
        esc += nesc > 0
        ncs.append(len(at)); trees += (len(at) == 1 and nesc == 0)
        for c in at:
            mins.append(c[0]); maxs.append(c[2]); lens.append(c[1])
    print(f"   p={p:<6} N={NN:<8} seeds {seeds}: trees {trees} ({100*trees/seeds:.1f}%), mean attracting cycles {statistics.mean(ncs):.2f},"
          f" max {max(ncs)}; cycle minima: median {statistics.median(mins):.0f}, 90% {sorted(mins)[int(.9*len(mins))]}, max {max(mins)};"
          f" largest cycle element {max(maxs)}; longest cycle {max(lens)}; maps with escapes {esc}")
print(f"   ({time.time()-t0:.1f}s)")
# long cycles of random members: recompute the up-count a with the same hash, in Python
M64 = (1 << 64) - 1
def smix(x):
    x = (x + 0x9E3779B97F4A7C15) & M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
    return x ^ (x >> 31)
longc = set()
for p, seeds, NN in [(0.5, 300, 10**5), (0.1, 300, 10**5), (0.01, 300, 10**5)]:
    thr = int(p * 2.0**64)
    for s_ in range(1, seeds + 1):
        nesc, cyc = parse_cycles(run("random", s_, p, 0, NN, CAP)[0])
        for mn, ln, mx, bs in cyc:
            if ln >= 20:
                e = lambda i: 1 if smix(((s_ * 0x100000001B3) & M64) ^ smix(i)) < thr else 0
                v, a = mn, 0
                for _ in range(ln):
                    i = (v + 1) // 2; up = (v & 1) ^ e(i); a += up; v = v + i if up else v - i
                assert v == mn
                longc.add((a, ln))
print("   cycles of length >= 20 in these random members, as (up-steps a, length K):", sorted(longc))
print("   K/a:", sorted(set(round(K / a, 5) for a, K in longc)), " (log2 3 = 1.58496)")
print(f"   max |K - a log2 3| over these cycles: {max(abs(K - a * math.log2(3)) for a, K in longc):.3f}"
      "  (2^K and 3^a agree within that many bits: every long cycle sits at a good approximation of log2 3)")

# ------------------------------------------------------------------------------------------ (F)
print("\n(F) Single flips: change exactly one pair of T (offset 0) or of U = 3n-1 (offset 1)")
for off in (0, 1):
    out = run("single", off, 10**7, 10**6)
    frag = [l for l in out if l.startswith("FRAG")]
    summ = [l for l in out if l.startswith("SINGLE")]
    idx = [int(l.split()[2]) for l in frag]
    print(f"   offset {off}: {summ[-1]}")
    print(f"     fragile pairs i (a new cycle appears): {idx}")
    ann = []
    for l in frag:
        f = l.split(); i = int(f[2]); mn = int(f[6])
        pa = (2 * i - 1, 2 * i) if off == 0 else (2 * i, 2 * i + 1)
        def Ff(n):
            j = (n + 1) // 2 if off == 0 else n // 2
            up = (n & 1) ^ (1 if j == i else 0)
            return n + j if up else n - j
        v, a, K = mn, 0, 0
        while True:
            w = Ff(v); a += w > v; K += 1; v = w
            if v == mn: break
        ann.append((mn, K, a, int(f[10])))
    print(f"     new cycles (min, length K, up-steps a, max): {ann}")
    print(f"     distinct (a, K): {sorted(set((a, K) for _, K, a, _ in ann))}  (3^a ~ 2^K: a/K near log2/log3 = 0.6309)")
print("   Mechanism: after flipping pair i the new cycle is a path 3i ~> 2i (or i-1 ~> 2i-1) of the old map, i.e.")
print("   i(2^(K+1) - 3^(a+1)) = B_w: the same Diophantine gate as a cycle, so fragile pairs thin out like cycles do")
print("   (they may recur near later approximations of log2 3, as 84/53 does for 3n-1 at i = 12029).")

# ------------------------------------------------------------------------------------------ (G)
print("\n(G) The provable rung: 'landing => down' (every up-move lands on a down-mover)")
print("    PROVED: then n -> ceil(3n/2) -> floor(ceil(3n/2)/2) < n for n >= 2, so every orbit reaches the root: a tree.")
print("    Constraints (derived): eps_3k = 1 - eps_2k;  eps_(2k+1) = 0 => eps_(3k+1) = 0;  eps_(2k+1) = 1 => eps_(3k+2) = 1.")
print("    Collatz violates it at every even i (3i-1 odd: the rising runs).  No periodic member satisfies it (section D).")
for l in run("landing", 2**24, 10**6):
    print("   ", l)
for l in run("landingvc", 10**7):
    print("   ", l)
try:
    from ortools.sat.python import cp_model
    def solve(X):
        m = cp_model.CpModel(); top = (3 * X) // 2 + 4
        e = [m.NewBoolVar("") for _ in range(top + 1)]
        m.Add(e[1] == 0)
        for i in range(1, X + 1):
            k = i // 2
            if i % 2 == 0:
                m.Add(e[2 * k] + e[3 * k] == 1)
            else:
                if 3 * k + 1 != i: m.AddImplication(e[3 * k + 1], e[i])
                m.AddImplication(e[i], e[3 * k + 2])
        m.Minimize(sum(e[1:X + 1]))
        s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 120
        s.parameters.max_memory_in_mb = 600
        st = s.Solve(m)
        return s.StatusName(st), int(s.ObjectiveValue()), int(s.BestObjectiveBound())
    for X in (3000, 30000, 100000):
        st, ob, bd = solve(X)
        print(f"    CP-SAT exact window optimum, sources i <= {X}: status {st}, min flips {ob} (bound {bd}), density {ob/X:.5f}")
except Exception as ex:  # pragma: no cover
    print("    (ortools unavailable:", ex, ")")
print("    => every member of the landing family differs from Collatz on >= ~29% of the pairs in [1,1e5] (exact window")
print("       optimum), against 1/4 from the vertex-cover relaxation and 1/3 for the greedy member.")

# ------------------------------------------------------------------------------------------ (H)
print("\n(H) A density-zero modification with a divergent orbit (PROVED construction)")
chain = [3]
while chain[-1] < 10**30:
    n = chain[-1]; chain.append(n + (n + 1) // 2)
flipped = [((n + 1) // 2) for n in chain if n % 2 == 0]      # Collatz moves even n down; flip their pairs
pairs = [(n + 1) // 2 for n in chain]
distinct = len(set(pairs)) == len(pairs)
eps = set(flipped)
def Fchain(n):
    i = (n + 1) // 2; up = (n & 1) ^ (1 if i in eps else 0)
    return n + i if up else n - i
v, ok = 3, True
for t in range(len(chain) - 1):
    if Fchain(v) != chain[t + 1]: ok = False; break
    v = Fchain(v)
print(f"   chain c_0 = 3, c_(t+1) = c_t + ceil(c_t/2): {chain[:12]} ...  ({len(chain)} terms below 1e30)")
print(f"   pairs of the chain are distinct: {distinct}; flipped pairs (even chain members): {len(flipped)} of {len(chain)};")
print(f"   the modified map follows the chain upward (verified to 1e30): {ok}")
print("   The flip set has at most log_(3/2) X + 1 elements below X: density 0, yet the orbit of 3 diverges.")
print("   Every such map still preserves all pair sums and uses every length once up and once down.")
# how visible is such a planted chain to orbit statistics?  fraction of n <= NN whose modified orbit reaches the chain
def captured(c0, NN):
    ch = [c0]
    while ch[-1] < 10**40:
        ch.append(ch[-1] + (ch[-1] + 1) // 2)
    chs = set(ch); fl = {(x + 1) // 2 for x in ch if x % 2 == 0}
    st = bytearray(NN + 1); div = 0
    for n0 in range(1, NN + 1):
        if not st[n0]:
            path, v, res = [], n0, 0
            while True:
                if v in chs: res = 2; break
                if v <= NN and st[v]: res = st[v]; break
                if v == 1: res = 1; break
                path.append(v)
                i = (v + 1) // 2; up = (v & 1) ^ (1 if i in fl else 0)
                v = v + i if up else v - i
            for p in path:
                if p <= NN: st[p] = res
        div += st[n0] == 2
    return div / NN
print("   orbit statistics do see a planted chain, with weight falling roughly like 1/c_0: fraction of n <= 1e5 captured:")
print("     " + ", ".join(f"c_0 = {c0}: {captured(c0, 10**5):.5f}" for c0 in (3, 27, 1001, 100003)))
print("   (the chain from 3 runs through 8 = 2*4, the pair {7,8}, which carries almost every Collatz orbit)")
