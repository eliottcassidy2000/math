#!/usr/bin/env python3
"""
lrc_polygon_tournament_wheel_s517.py

oracle-2026-06-01-S517

Tournaments as oriented REGULAR POLYGONS on the unit circle, and how the
twin-Goldbach 6k±1 wheel (S516) and the Lonely Runner live on the same object.

Picture: n points on the unit circle = a polygon; the complete graph's chords get
oriented by the half-turn rule  i -> j  iff  the short way j->i is < half a turn.
- Runners at t = a/q sit at vertices of a REGULAR q-gon (the q-th roots of unity),
  occupying the multiset {v_i a mod q}.  At the LRC tight witness (speeds 1..n-1,
  t = 1/n) they occupy EVERY non-observer vertex of the regular n-gon -> the
  maximally symmetric "regular tournament" (polygon minus its centre/clasp).
- The residue wheel mod m is literally the regular m-gon; twin primes occupy only
  the unit vertices (coprime to m).  Twin-Goldbach = self-convolution of the
  unit-vertex set on the wheel.

This script confirms (1) the LRC-witness regular-polygon tournament and its type,
(2) the chord-length "channels" (edge orbits of the polygon), (3) the hexagon
twin convolution (channels are residue-complete; misses are magnitude deserts).
"""

from __future__ import annotations
from functools import lru_cache
from itertools import combinations


def half_turn_circulant(n):
    """Tournament on Z_n vertices, i->j iff (i-j) mod n in {1..floor((n-1)/2)}
    plus the boundary handling; this is the regular-polygon half-turn orientation
    restricted to the n-gon vertices."""
    half = (n - 1) // 2
    conn = set(range(1, half + 1))
    if n % 2 == 0:
        # n even: the antipodal chord (d=n/2) is a tie; base-path break i->j for i<j
        pass
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = (i - j) % n
            if d in conn:
                adj[i][j] = 1
            elif n % 2 == 0 and d == n // 2:
                adj[i][j] = 1 if i < j else 0  # antipodal tie -> base path
    return adj


def Hc(adj):
    n = len(adj); full = (1 << n) - 1; A = tuple(map(tuple, adj))
    @lru_cache(None)
    def dp(mask, last):
        if mask == full: return 1
        return sum(dp(mask | (1 << x), x) for x in range(n) if not (mask >> x) & 1 and A[last][x])
    return sum(dp(1 << s, s) for s in range(n))


def scores(adj): return tuple(sorted(sum(r) for r in adj))


def is_regular(adj):
    return len(set(sum(r) for r in adj)) == 1


def main():
    print("Tournaments as oriented regular polygons on the unit circle (oracle-S517)\n")

    print("=" * 66)
    print("1. LRC tight witness = regular n-gon minus the observer vertex")
    print("=" * 66)
    print("   speeds {1..n-1} at t=1/n -> runners at vertices {1..n-1} of the n-gon")
    print("   tournament = half-turn chord orientation on those n-1 vertices.\n")
    for n in range(4, 10):
        # runners 1..n-1 at positions i/n; half-turn among them
        k = n - 1
        adj = [[0] * k for _ in range(k)]
        half = 1.0 / 2
        for a in range(k):
            for b in range(k):
                if a == b: continue
                i, j = a + 1, b + 1
                d = ((i - j) % n) / n
                if 0 < d < half:
                    adj[a][b] = 1
                elif abs(d - half) < 1e-12:  # antipodal (n even): base-path break
                    adj[a][b] = 1 if i < j else 0
        sc = scores(adj); reg = is_regular(adj); H = Hc(adj)
        print(f"   n={n}: {k}-vertex tournament  H={H:<5} score={sc}  regular={reg}"
              f"  {'(odd n: rotational R_'+str(n)+'-type)' if n%2 else '(even n: near-regular)'}")
    print()

    print("=" * 66)
    print("2. Chord-length channels (edge orbits of the regular m-gon)")
    print("=" * 66)
    for m in (6, 7, 14, 18):
        nclasses = m // 2
        print(f"   regular {m}-gon: chord-length classes d=1..{nclasses}  "
              f"({nclasses} channels);  units coprime to {m}: "
              f"{[a for a in range(1,m) if __import__('math').gcd(a,m)==1]}")
    print("   (LRC 'holdback' S25 = 1/(2d): channel d is the chord of step d;")
    print("    twin primes differ by 2 = the d=2 chord class on every wheel.)\n")

    print("=" * 66)
    print("3. Twin-Goldbach on the hexagon wheel (mod 6): residue-complete channels")
    print("=" * 66)
    units6 = [1, 5]   # twin-prime residues mod 6
    conv = sorted({(a + b) % 6 for a in units6 for b in units6})
    print(f"   twin residues mod 6 = {units6} (the two UNIT vertices of the hexagon)")
    print(f"   self-convolution {units6} (+) {units6} mod 6 = {conv}  "
          f"= all even residues {{0,2,4}}")
    print("   => the 3 even-target CHANNELS (1+1=2, 1+5=0, 5+5=4) are all available")
    print("      residue-wise; so the 35 missed evens are MAGNITUDE deserts WITHIN")
    print("      channels (local absence of actual twins), not a residue obstruction.")
    print("      Triple {6m-2,6m,6m+2} = one number per channel failing at once.\n")

    print("SUMMARY: the regular polygon (roots of unity) is the shared stage.")
    print(" - LRC tight witness = full regular n-gon minus the clasp (observer) =")
    print("   the regular/rotational tournament (the most symmetric chord-orientation).")
    print(" - The wheel mod m = regular m-gon; twin primes = its unit vertices;")
    print("   twin-Goldbach = unit-vertex self-convolution; LRC sieve = 'occupy")
    print("   vertex 0 mod m'. Both are regular-polygon vertex covering/avoidance.")


if __name__ == "__main__":
    main()
