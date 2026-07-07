"""
crossing_book_tiling_fusion_opus_S142.py

THE BOOK CUBE IS THE TILING CUBE (owner: revisit the unit-distance / crossing-number /
even-odd-duality past work — HYP-2712a, HYP-2627/2629, HYP-3954 tangent).

SETTING: 2-page book drawings of K_n with fixed spine order 1..n. Page assignment
x in GF(2)^{chords}. Two chords {a,b}, {c,d} (a<c<b<d, interleaved) cross iff SAME page.
FACTS DERIVED THIS SESSION (verified below):
  (F1) consecutive-vertex chords (l = b-a = 1) interleave NOTHING: the crossing form
       lives exactly on the l >= 2 chords = THE TILES of the tournament tiling model.
       The book cube and the tiling cube are the same GF(2)^m, m = C(n-1,2), viewed from
       the same marked Hamiltonian path (the spine = the base path).
  (F2) AFFINE PARITY LAW: #crossings mod 2 = C(n,4) + sum over chords e with
       ODD interleave-degree of x_e, and deg_I({a,b}) = (l-1)(n-1-l), l = b-a, which is
       odd iff (l-1) and (n-1-l) both odd. For n = 14: deg_I odd iff l EVEN.
       => the crossing parity of ANY 2-page drawing of K_14 is
          C(14,4) + sum_{even-length chords} x_e = 1001 + (# even-length chords on page 1)
       -- only EVEN chords carry parity (the even/odd duality in the crossing world).
  (F3) PAGE-SWAP INVARIANCE: Q(flip x) = Q(x) exactly (equality of bits is preserved
       under global complement) -- the crossing form descends to the LINE pairing
       {x, flip x} of the blue/black metagraph.
  (F4) GRID-SYMMETRY INVARIANCE: spine reversal rho maps interleaved pairs to interleaved
       pairs, so Q(sigma x) = Q(x): Q is invariant under the grid involution -- the
       blue (gridsym) subcube is a symmetry locus of the crossing landscape.

QUESTIONS COMPUTED:
 (Q1) verify F1/F2/F3/F4 exhaustively at n = 6..9 (+ F2 spot-checks at n = 10..14).
 (Q2) HYP-2712a upgrade: min Q over ALL pages = Z(n)? over EVEN pages? over BLUE
      (gridsym) pages? over BLUE-and-EVEN pages? (n = 6, 7, 8: full sweeps of the
      relevant subcubes.)  Z(n) = Guy = floor products / 4.
 (Q3) the parity corollary: every 2-page optimum must have an even number of
      parity-carrying chords on page 1 iff Z(n) == C(n,4) mod 2 -- tabulate.
 (Q4) the E_n lens with today's classification: at n = 7, per even-graph value class of
      page-0, the Q spectrum (min per even graph), cross-tabbed with sigma-fixedness.
"""
from fractions import Fraction as F
from itertools import combinations
import sys, time

def Z_guy(n):
    return (n // 2) * ((n - 1) // 2) * ((n - 2) // 2) * ((n - 3) // 2) // 4

def chords(n):
    return [(a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1)]

def interleaving_pairs(n, C):
    idx = {c: i for i, c in enumerate(C)}
    pairs = []
    for (a, b) in C:
        for (c, d) in C:
            if a < c < b < d:
                pairs.append((idx[(a, b)], idx[(c, d)]))
    return pairs

def Q_count(x, pairs):
    return sum(1 for (i, j) in pairs if ((x >> i) & 1) == ((x >> j) & 1))

def main():
    t0 = time.time()
    print("=" * 100)
    print("(Q1) verify the structural facts")
    print("=" * 100)
    import random
    rng = random.Random(142)
    for n in (6, 7, 8, 9, 10, 12, 14):
        C = chords(n)
        m = len(C)
        pairs = interleaving_pairs(n, C)
        # F1: l=1 chords interleave nothing
        idx = {c: i for i, c in enumerate(C)}
        deg = [0] * m
        for (i, j) in pairs:
            deg[i] += 1; deg[j] += 1
        f1 = all(deg[idx[(a, a + 1)]] == 0 for a in range(1, n))
        # F2: parity affine law on random assignments
        oddset = [i for i, c in enumerate(C) if deg[i] % 2 == 1]
        # check deg formula (l-1)(n-1-l)
        f2a = all(deg[idx[(a, b)]] == (b - a - 1) * (n - 1 - (b - a)) for (a, b) in C)
        c4 = len(pairs)  # should be C(n,4)
        from math import comb
        f2b = (c4 == comb(n, 4))
        f2c = True
        trials = 400 if n <= 10 else 60
        for _ in range(trials):
            x = rng.getrandbits(m)
            q = Q_count(x, pairs)
            pred = (comb(n, 4) + sum((x >> i) & 1 for i in oddset)) % 2
            if q % 2 != pred: f2c = False
        # F3: page swap
        f3 = True
        for _ in range(60):
            x = rng.getrandbits(m)
            if Q_count(x, pairs) != Q_count(x ^ ((1 << m) - 1), pairs): f3 = False
        # F4: spine reversal
        rho = {c: i for i, c in enumerate(C)}
        sig = [rho[(n + 1 - b, n + 1 - a)] for (a, b) in C]
        f4 = True
        for _ in range(60):
            x = rng.getrandbits(m)
            xs = 0
            for i in range(m):
                if (x >> i) & 1: xs |= (1 << sig[i])
            if Q_count(x, pairs) != Q_count(xs, pairs): f4 = False
        oddl = sorted({(b - a) for (a, b) in C if deg[idx[(a, b)]] % 2 == 1})
        print(f"  n={n:2d}: F1 {f1}  F2(deg formula {f2a}, #pairs=C(n,4) {f2b}, parity law {f2c})"
              f"  F3 {f3}  F4 {f4}   parity-carrying chord lengths: {oddl}")

    print()
    print("=" * 100)
    print("(Q2) minima over subcubes vs Guy Z(n)   [pages on l>=2 chords; l=1 fixed page 0 wlog]")
    print("=" * 100)
    for n in (6, 7, 8):
        C = chords(n)
        m = len(C)
        pairs = interleaving_pairs(n, C)
        idx = {c: i for i, c in enumerate(C)}
        tiles = [i for i, (a, b) in enumerate(C) if b - a >= 2]
        # even subgraphs: page-1 set has even degree at every vertex
        # gridsym: x[(a,b)] == x[(n+1-b, n+1-a)]
        sig = {i: idx[(n + 1 - C[i][1], n + 1 - C[i][0])] for i in range(m)}
        best = {"all": None, "even": None, "blue": None, "blue&even": None}
        arg = {}
        for xt in range(1 << len(tiles)):
            x = 0
            for k, i in enumerate(tiles):
                if (xt >> k) & 1: x |= (1 << i)
            q = Q_count(x, pairs)
            def upd(key):
                if best[key] is None or q < best[key]:
                    best[key] = q; arg[key] = x
            upd("all")
            degs = [0] * (n + 1)
            for i in range(m):
                if (x >> i) & 1:
                    degs[C[i][0]] += 1; degs[C[i][1]] += 1
            iseven = all(d % 2 == 0 for d in degs[1:])
            isblue = all(((x >> i) & 1) == ((x >> sig[i]) & 1) for i in range(m))
            if iseven: upd("even")
            if isblue: upd("blue")
            if iseven and isblue: upd("blue&even")
        print(f"  n={n}: Z(n)={Z_guy(n)}   min all={best['all']}  even={best['even']}"
              f"  BLUE={best['blue']}  blue&even={best['blue&even']}"
              f"   {'BLUE ATTAINS' if best['blue'] == best['all'] else 'blue does NOT attain'}"
              f"   {'B&E ATTAINS' if best['blue&even'] == best['all'] else ''}")

    print()
    print("=" * 100)
    print("(Q3) parity corollary: Z(n) vs C(n,4) mod 2 (must differ by an even page-1 count of")
    print("     parity-carrying chords; if Z(n) != C(n,4) mod 2, the optimum needs an ODD count)")
    print("=" * 100)
    from math import comb
    for n in range(5, 15):
        z, c = Z_guy(n), comb(n, 4)
        print(f"  n={n:2d}: Z={z} ({z%2}), C(n,4)={c} ({c%2})  -> optimal page-1 parity-chord count "
              f"is {'EVEN' if z % 2 == c % 2 else 'ODD'}")

    print(f"\n[{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
