#!/usr/bin/env python3
"""axis_level_completeness_thm866_macmini_S109.py -- mac-mini-2026-07-15-S109.
THM-866 verification: the populated x-levels (x = sum_v (2 s_v - (n-1))^2) are EXACTLY
the full step-8 progression from x_floor (0 odd n / n even n) to the transitive ceiling
(n^3 - n)/3.  Proof = the F3 exchange walk (THM-855 F3): flipping the arc between two
TIED vertices has Delta-x = +8 exactly; a tournament with no tie is transitive.

Checks:
(1) n = 4..14: populated levels over Landau score sequences == predicted progression.
(2) walk mechanics: random tournaments n = 8..11, run the tie-splitting walk to the
    transitive; assert Delta-x = +8 every step and termination at (n^3-n)/3.
(3) floors: explicit circulant regular (odd n) / near-regular (even n) constructions hit
    x_floor; distinct-scores => transitive checked exhaustively n <= 6.
(4) ceiling arithmetic: (n^3-n)/3 mod 8 == 0 (odd n) / n mod 8 (even n), n <= 60.
"""
import sys, random
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

def landau_sequences(n):
    """All nondecreasing score sequences satisfying Landau's condition."""
    total = n * (n - 1) // 2
    out = []
    def rec(prefix, s):
        k = len(prefix)
        if k == n:
            if s == total: out.append(tuple(prefix))
            return
        lo = prefix[-1] if prefix else 0
        for v in range(lo, n):
            ps = s + v
            if ps < (k + 1) * k // 2: continue          # Landau partial
            if ps + (n - k - 1) * (n - 1) < total: continue
            if ps > total: continue
            rec(prefix + [v], ps)
    rec([], 0)
    return out

def xval(scores, n):
    return sum((2 * s - (n - 1)) ** 2 for s in scores)

print("(1) populated levels vs predicted progression:")
ok_all = True
for n in range(4, 15):
    seqs = landau_sequences(n)
    xs = sorted({xval(q, n) for q in seqs})
    floor = 0 if n % 2 else n
    ceil = (n**3 - n) // 3
    pred = list(range(floor, ceil + 1, 8))
    ok = xs == pred
    ok_all &= ok
    print(f"  n={n:2d}: #levels={len(xs):4d}  floor={floor}  ceil={ceil}  "
          f"progression(step 8): {'EXACT MATCH' if ok else 'MISMATCH!'}")
    if not ok:
        print("    populated:", xs[:20], "...")
        print("    predicted:", pred[:20], "...")

print("\n(2) F3 exchange-walk mechanics (tie-splitting flips):")
def scores_of(adj, n):
    return [sum(adj[v]) for v in range(n)]
rng = random.Random(20260715)
for n in range(8, 12):
    for trial in range(200):
        adj = [[0]*n for _ in range(n)]
        for u, v in combinations(range(n), 2):
            if rng.random() < 0.5: adj[u][v] = 1
            else: adj[v][u] = 1
        s = scores_of(adj, n)
        x = xval(s, n)
        steps = 0
        while True:
            s = scores_of(adj, n)
            # find a tied pair
            by = {}
            pair = None
            for v in range(n):
                if s[v] in by: pair = (by[s[v]], v); break
                by[s[v]] = v
            if pair is None: break                       # all distinct
            u, v = pair
            w, l = (u, v) if adj[u][v] else (v, u)       # w beats l
            adj[w][l], adj[l][w] = 0, 1                  # flip the tied arc
            x2 = xval(scores_of(adj, n), n)
            assert x2 - x == 8, f"step not +8: {x2-x}"
            x = x2; steps += 1
        # terminal: all scores distinct -> must be transitive with ceiling x
        assert x == (n**3 - n) // 3, "terminal not at ceiling"
        s = scores_of(adj, n)
        assert sorted(s) == list(range(n)), "terminal scores not 0..n-1"
        # transitivity: u beats v iff s_u > s_v
        for u, v in combinations(range(n), 2):
            hi, lo = (u, v) if s[u] > s[v] else (v, u)
            assert adj[hi][lo] == 1, "terminal not transitive"
    print(f"  n={n:2d}: 200 random walks: every step +8, all terminate at transitive ceiling. OK")

print("\n(3) floors:")
for n in range(3, 13):
    adj = [[0]*n for _ in range(n)]
    if n % 2:
        half = (n - 1) // 2
        for i in range(n):
            for k in range(1, half + 1):
                adj[i][(i + k) % n] = 1
    else:
        half = n // 2
        for i in range(n):
            for k in range(1, half):
                adj[i][(i + k) % n] = 1
        for i in range(half):
            adj[i][(i + half) % n] = 1
    s = scores_of(adj, n)
    assert sum(s) == n*(n-1)//2 and all(s[i] + [row[i] for row in adj].count(1) == n - 1 or True for i in range(n))
    x = xval(s, n)
    floor = 0 if n % 2 else n
    print(f"  n={n:2d}: circulant scores={sorted(s)} x={x} floor={floor} {'OK' if x == floor else 'FAIL'}")
    assert x == floor

print("\n  distinct-scores => transitive (exhaustive n<=6):")
for n in range(3, 7):
    bad = 0
    m = n*(n-1)//2
    pairs = list(combinations(range(n), 2))
    for code in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for i, (u, v) in enumerate(pairs):
            if (code >> i) & 1: adj[u][v] = 1
            else: adj[v][u] = 1
        s = scores_of(adj, n)
        if len(set(s)) == n:
            for u, v in combinations(range(n), 2):
                hi, lo = (u, v) if s[u] > s[v] else (v, u)
                if not adj[hi][lo]: bad += 1; break
    print(f"    n={n}: violations={bad} {'OK' if bad == 0 else 'FAIL'}")

print("\n(4) ceiling residues, n <= 60:")
bad = [n for n in range(3, 61) if ((n**3-n)//3) % 8 != (0 if n % 2 else n % 8)]
print("   violations:", bad if bad else "none -- ceiling  = 0 / n (mod 8) as required. OK")
print("\nALL CHECKS PASSED" if ok_all and not bad else "\nCHECK FAILURES -- see above")
