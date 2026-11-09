"""Probe for an involution proof of the blue-line H-mod-4 law: for each sigma-orbit of tiles (fixed anti-diagonal tiles
and transposed pairs), flipping the orbit keeps a blue tiling blue.  Record, over all blue tilings t, whether the flip
changes H mod 4 (delta = 1) always, never, or sometimes.
usage: python blue_flip_probe.py n
"""
import sys
import itertools


def setup(n):
    tiles = [(i, j) for i in range(n) for j in range(i + 2, n)]
    idx = {t: k for k, t in enumerate(tiles)}
    sg = [idx[(n - 1 - j, n - 1 - i)] for (i, j) in tiles]
    orbits, seen = [], set()
    for k in range(len(tiles)):
        if k in seen:
            continue
        o = {k, sg[k]}
        seen |= o
        orbits.append(sorted(o))
    return tiles, orbits


def H(n, tiles, bits):
    out = [0] * n
    for i in range(n - 1):
        out[i] |= 1 << (i + 1)
    for k, (i, j) in enumerate(tiles):
        if bits >> k & 1:
            out[j] |= 1 << i
        else:
            out[i] |= 1 << j
    F = 1 << n
    f = [[0] * n for _ in range(F)]
    for v in range(n):
        f[1 << v][v] = 1
    for S in range(1, F):
        for v in range(n):
            c = f[S][v]
            if not c:
                continue
            o = out[v] & ~S
            while o:
                u = (o & -o).bit_length() - 1
                o &= o - 1
                f[S | (1 << u)][u] += c
    return sum(f[F - 1])


def main():
    n = int(sys.argv[1])
    tiles, orbits = setup(n)
    m = len(tiles)
    allones = (1 << m) - 1
    masks = [sum(1 << k for k in o) for o in orbits]
    blue = []
    for code in range(1 << len(orbits)):
        t = 0
        for q in range(len(orbits)):
            if code >> q & 1:
                t |= masks[q]
        blue.append(t)
    Hval = {t: H(n, tiles, t) for t in blue}
    print(f"n = {n}: {len(blue)} blue tilings, {len(orbits)} sigma-orbits "
          f"({sum(1 for o in orbits if len(o) == 1)} fixed tiles)")
    for q, o in enumerate(orbits):
        deltas = set()
        for t in blue:
            u = t ^ masks[q]
            deltas.add(((Hval[u] - Hval[t]) // 2) % 2)
        kind = "fixed" if len(o) == 1 else "pair"
        print(f"  orbit {q:2d} {kind:5s} tiles {[tiles[k] for k in o]}: delta(H mod 4) in {sorted(deltas)}")
    # joint distribution over blue lines
    cnt = {}
    for t in blue:
        if t & (1 << (m - 1)):
            continue
        a, b = Hval[t] % 4, Hval[t ^ allones] % 4
        key = tuple(sorted((a, b)))
        cnt[key] = cnt.get(key, 0) + 1
    print("  blue lines by (H mod 4) pair:", dict(sorted(cnt.items())))


if __name__ == "__main__":
    main()
