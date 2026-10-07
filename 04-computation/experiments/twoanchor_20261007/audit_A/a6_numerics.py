#!/usr/bin/env python3
"""Audit A, items 6/7: numerical claims of THM-4601, own code.

Pair chain (THM-4581 table, own derivation) in exact integer form: state (k, E), e = E / 3^max(0,-k).
  sigma = E mod 2, beta = parity bit of the base orbit v.
  (0,0): E/2                         (0,1): k>=0: (3E+1-3^k)/2      k<0: (3E+3^|k|-1)/2
  (1,0): k+1, k>=0: (3E+1)/2   k<0: (E+3^(|k|-1))/2
  (1,1): k-1, k>=1: (E-3^(k-1))/2   k<=0: (3E-1)/2
Absorbed <=> (k, E) = (0, 0).

A  exact BFS of absorbing bit patterns from the universal state (3, -26): counts by first-absorption length,
   unconditioned and conditioned on the post-run type (r = 1: first bits 1,1;  r = 2: first bits 1,0,0).
B  head-grammar completeness: exact conditional BFS mass by post-run depth 23 vs the session's union of
   (debt 3, i = 1) head classes with sum(u) <= 22  (0.013902 for r = 1, 0.029463 for r = 2).
C  Monte Carlo P(absorbed by s) from (3, -26): fair bits, and conditioned on r = 1 / r = 2.
D  re-anchoring at -1 after the two-run: the all-ones post-run pattern, on the chain and on an actual integer source.
E  small J-table (reset child D = 1 vs D = 3) on actual sources, 400-bit t, post-run budget 400.
"""
import random, sys, math
from collections import defaultdict, Counter
from a2_barrier import make_source, T

def step(k, E, beta):
    sigma = E & 1
    if sigma == 0 and beta == 0:
        return k, E >> 1
    if sigma == 0 and beta == 1:
        if k >= 0:
            return k, (3*E + 1 - 3**k) >> 1
        return k, (3*E + 3**(-k) - 1) >> 1
    if sigma == 1 and beta == 0:
        if k >= 0:
            return k + 1, (3*E + 1) >> 1
        return k + 1, (E + 3**(-k - 1)) >> 1
    if k >= 1:
        return k - 1, (E - 3**(k - 1)) >> 1
    return k - 1, (3*E - 1) >> 1

def selftest(n=3000):
    """the integer chain against actual orbits u = 3^k v + e"""
    rnd = random.Random(1)
    from fractions import Fraction as Fr
    for _ in range(n):
        v = rnd.getrandbits(200) | 1
        k = 3; E = -26
        u = 27*v - 26
        for s in range(150):
            beta = v & 1
            k, E = step(k, E, beta)
            u, v = T(u), T(v)
            e = Fr(E, 3**max(0, -k))
            assert Fr(u) == Fr(3)**k * v + e, "chain mismatch"
            if k == 0 and E == 0:
                assert u == v
                break
    print(f"   self-test: integer chain = actual orbit relation on {n} random 200-bit pairs (150 steps)")

def bfs(prefix=(), depth=23, maxstates=6_000_000):
    """counts[L] = number of length-L bit strings (with given prefix) whose first absorption is at time L."""
    cur = {(3, -26): 1}
    counts = Counter()
    L = 0
    for L in range(1, depth + 1):
        nxt = defaultdict(int)
        bits = (prefix[L - 1],) if L <= len(prefix) else (0, 1)
        for (k, E), c in cur.items():
            for b in bits:
                k2, E2 = step(k, E, b)
                if k2 == 0 and E2 == 0:
                    counts[L] += c
                else:
                    nxt[(k2, E2)] += c
        cur = nxt
        if len(cur) > maxstates:
            return counts, L, len(cur)
    return counts, L, len(cur)

def part_A_B():
    counts, L, ns = bfs((), 20)
    seq = [counts[l] for l in range(1, 21)]
    print(f"A  unconditioned absorbing-pattern counts, lengths 1..20: {seq}")
    print(f"   lengths 9..20: {seq[8:]}  (session: 1, 1, 3, 7, 15, 30, 67, 147, 301, 658, 1357, 2783)")
    mass20 = sum(counts[l] / 2**l for l in counts)
    print(f"   unconditioned mass absorbed by 20: {mass20:.6f}; states alive at 20: {ns}")
    # which first patterns of length 9..11 ?
    res = {}
    for r, pre in ((1, (1, 1)), (2, (1, 0, 0)), ('letter2', (1, 0, 1)), ('even', (0,))):
        c, Lr, nsr = bfs(pre, 23)
        m = {d: sum(c[l] for l in c if l <= d) / 2**(d - len(pre)) for d in (12, 20, 23)}
        mass = lambda d: sum(c[l] / 2**(l - len(pre)) for l in c if l <= d)
        first = min(c) if c else None
        res[r] = (first, mass(20), mass(23), Lr, nsr)
        print(f"   prefix {pre} ({r}): first absorption at length {first}; conditional mass by 20: {mass(20):.6f}, "
              f"by 23: {mass(23):.6f} (BFS reached depth {Lr}, {nsr} live states)")
    print("B  head-grammar check (session: union of debt-3 heads with sum(u) <= 22 = absorption by post-run depth 23):")
    print(f"   r=1: exact conditional BFS mass by 23 = {res[1][2]:.6f}  vs session union 0.013902")
    print(f"   r=2: exact conditional BFS mass by 23 = {res[2][2]:.6f}  vs session union 0.029463")
    return res

def mc(nrun, smax, prefix, rnd, checkpoints):
    hits = Counter()
    for _ in range(nrun):
        k, E = 3, -26
        absorbed_at = None
        for s in range(smax):
            b = prefix[s] if s < len(prefix) else rnd.getrandbits(1)
            k, E = step(k, E, b)
            if k == 0 and E == 0:
                absorbed_at = s + 1; break
        if absorbed_at is not None:
            for c in checkpoints:
                if absorbed_at <= c:
                    hits[c] += 1
    return {c: hits[c]/nrun for c in checkpoints}

def part_C(nrun=20000, smax=1000):
    rnd = random.Random(20261007)
    cps = (30, 100, 200, 300, 400, 1000)
    for name, pre in (("fair bits", ()), ("r=1 (1,1)", (1, 1)), ("r=2 (1,0,0)", (1, 0, 0))):
        p = mc(nrun, smax, pre, rnd, cps)
        se = {c: math.sqrt(p[c]*(1-p[c])/nrun) for c in cps}
        print(f"C  {name:12s} P(absorbed by s), {nrun} runs: " + ", ".join(f"s={c}: {p[c]:.3f}±{se[c]:.3f}" for c in cps))

def part_C_long(nrun=4000, smax=4000):
    rnd = random.Random(77)
    p = mc(nrun, smax, (), rnd, (1000, 4000))
    print(f"C' fair bits, {nrun} runs: P(by 1000) = {p[1000]:.3f}, P(by 4000) = {p[4000]:.3f}")

def part_D():
    # chain: all-ones bits after the universal state
    k, E = 3, -26
    traj = []
    for s in range(30):
        k, E = step(k, E, 1)
        traj.append((s + 1, k, E))
    anchored = [(s, k, E) for (s, k, E) in traj if k < 0 and E == 1 - 3**(-k)]
    print(f"D  all-ones post-run bits from (3,-26): states at s=12..14: {traj[11:14]}; "
          f"first -1-anchored state (E = 1 - 3^|k| scaled, i.e. u+1 = 3^k (v+1)): {anchored[:1]}")
    # actual source: r = 1 type, Y_end = -1 mod 2^40  (w' = -1 mod 2^39)
    rnd = random.Random(3)
    for (K, J) in ((10, 3), (50, 7), (400, 20)):
        r = 1; j = 2*J + r; M = j + 400
        wp = ((1 << 39) * (rnd.getrandbits(300) * 2 + 1) - 1)        # w' = -1 mod 2^39, odd
        ux = (pow(3, -(J - 3), 1 << (M - j)) * wp) % (1 << (M - j))
        t = (pow(3, -(K - 1), 1 << M) * (1 + (ux << (j - 1)))) % (1 << M)
        x = 2*3**(K - 1)*t - 1; y = (x + 1)//27 - 1
        u, v, kk = x, y, 3
        for _ in range(2*J):
            kk += (u & 1) - (v & 1); u, v = T(u), T(v)
        assert kk == 3 and u - 1 == 27*(v - 1) and (v + 1) % (1 << 40) == 0
        rel = []
        for s in range(1, 40):
            kk += (u & 1) - (v & 1); u, v = T(u), T(v)
            rel.append((s, kk, (u + 1)*3**(-kk) == v + 1 if kk < 0 else None))
        ok = [s for (s, kc, flag) in rel if flag]
        print(f"   actual source K={K}, J={J}, Y_end = -1 mod 2^40: debt after 12 post-run steps = {rel[11][1]}; "
              f"u+1 = 3^k (v+1) holds at post-run times {ok[0]}..{ok[-1]} (debt constant {set(rel[s-1][1] for s in ok)})")
    # reachability of -1-anchored states (k != 0) within 16 bits
    cur = {(3, -26): 1}
    found = Counter()
    for L in range(1, 17):
        nxt = defaultdict(int)
        for (k, E), c in cur.items():
            for b in (0, 1):
                k2, E2 = step(k, E, b)
                if k2 == 0 and E2 == 0:
                    continue
                if k2 != 0 and ((k2 > 0 and E2 == 3**k2 - 1) or (k2 < 0 and E2 == 1 - 3**(-k2))):
                    found[(L, k2)] += c
                nxt[(k2, E2)] += c
        cur = nxt
    print(f"   -1-anchored states (k != 0) reachable from (3,-26) within 16 bits, (length, k): count of patterns: "
          f"{sorted(found.items())[:10]}")

def part_E(N=150, S=400):
    rnd = random.Random(11)
    print(f"E  J-table (actual sources, 400-bit t, K in [12,40], {N} per row, post-run budget {S}): fraction certified")
    for J in (1, 3, 10, 16, 32):
        row = []
        for r in (1, 2):
            j = 2*J + r
            c1 = c3 = 0
            for _ in range(N):
                K = rnd.randint(12, 40)
                n, t, x = make_source(rnd, K, j, 400)
                for D in (1, 3):
                    y = (x + 1)//3**D - 1
                    u, v, k = x, y, D
                    hit = False
                    for s in range(2*J + S):
                        if u == v and k == 0:
                            hit = True; break
                        k += (u & 1) - (v & 1); u, v = T(u), T(v)
                    if D == 1: c1 += hit
                    else: c3 += hit
            row.append((c1/N, c3/N))
        print(f"   J={J:2d}: D=1 (reset child) r=1/r=2: {row[0][0]:.3f}/{row[1][0]:.3f}   D=3: {row[0][1]:.3f}/{row[1][1]:.3f}")

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "ABCDE"
    selftest()
    if "A" in which or "B" in which: part_A_B()
    if "C" in which: part_C(); part_C_long()
    if "D" in which: part_D()
    if "E" in which: part_E()
    print("DONE")
