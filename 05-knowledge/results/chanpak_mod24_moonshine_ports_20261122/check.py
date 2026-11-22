"""All 256 cuts and 36 physical arc seams at one nine-vertex source."""
from collections import Counter
from itertools import combinations
from math import gcd

N = 9
CODE = 0x7103D59DC
EDGES = list(combinations(range(N), 2))
FULL = (1 << N) - 1


def bits(mask):
    while mask:
        b = mask & -mask
        yield b.bit_length() - 1
        mask -= b


def endings(prev, start=None):
    dp = [[0] * N for _ in range(1 << N)]
    for t in range(N):
        dp[1 << t][t] = int(start is None or start == t)
    for subset in range(1, FULL + 1):
        if subset.bit_count() == 1:
            continue
        for t in bits(subset):
            rem = subset ^ (1 << t)
            dp[subset][t] = sum(dp[rem][k] for k in bits(rem & prev[t]))
    return dp


def values_at(code, cut, all_seams=True, cycles=False):
    incoming = [0] * N
    outgoing = [0] * N
    signs = []
    for j, (a, b) in enumerate(EDGES):
        forward = ((code >> j) ^ (cut >> a) ^ (cut >> b)) & 1
        u, v = (a, b) if forward else (b, a)
        incoming[v] |= 1 << u
        outgoing[u] |= 1 << v
        signs.append(1 if forward else -1)
    if cycles:
        dp = endings(incoming, start=0)
        return sum(dp[FULL][v] for v in bits(incoming[0]))
    left = endings(incoming)
    h = sum(left[FULL])
    if not all_seams:
        return h
    right = endings(outgoing)
    cover = [[0] * N for _ in range(N)]
    for subset in range(1, FULL):
        comp = FULL ^ subset
        for a in bits(subset):
            count = left[subset][a]
            if count:
                for b in bits(comp):
                    cover[a][b] += count * right[comp][b]
    return h, [s * (cover[b][a] - cover[a][b]) for s, (a, b) in zip(signs, EDGES)]


def transform(values, sign):
    result = values.copy()
    for i in range(N - 1):
        for s in range(len(result)):
            if s & (1 << i):
                result[s] += sign * result[s ^ (1 << i)]
    return result


def face_reconstruction(table, edge, modulus):
    root = next(i for i in range(N) if i not in edge)
    labels = [i for i in range(N) if i != root]
    def ordinary_cut(local):
        cut = sum(1 << labels[i] for i in bits(local))
        if cut & 1:
            cut ^= FULL
        return cut >> 1
    full = [table[ordinary_cut(s)] for s in range(1 << (N - 1))]
    marked = sum(1 << labels.index(i) for i in edge)
    def needed(s):
        if (s & marked) == marked:
            return False
        if modulus:
            return s.bit_count() <= 4
        return s.bit_count() <= 6 and (s.bit_count() < 6 or s & marked)
    samples = [v if needed(s) else 0 for s, v in enumerate(full)]
    coef = transform(samples, -1)
    recovery = transform([v if needed(s) else 0 for s, v in enumerate(coef)], 1)
    assert all(((a-b) % modulus == 0 if modulus else a == b)
               for a, b in zip(recovery, full))
    return sum(bool(needed(s)) for s in range(len(table)))


def main():
    h, d = zip(*(values_at(CODE, s << 1) for s in range(1 << (N - 1))))
    tables = list(zip(*d))
    for edge, table in zip(EDGES, tables):
        assert face_reconstruction(list(table), edge, 24) == 141
        assert face_reconstruction(list(table), edge, None) == 189
    print("All 36 seam tables recovered on 256 cuts, with 141 samples modulo 24 and 189 integer samples.")
    print("Orbit count residues modulo 24:", sorted(Counter(v % 24 for v in h).items()))
    moduli = [m for m in range(2, 514, 2) if len({v % m for v in h}) == m // 2]
    print("Even moduli with full odd support:", moduli)
    centers = Counter()
    for s in range(len(h)):
        for j, delta in enumerate(d[s]):
            if abs(delta) == 2:
                centers[h[s] + delta // 2] += 1
    print("Absolute-gap-two contacts:", sum(centers.values()), "distinct centers:", len(centers))
    print("Center support modulo 24:", sorted({c % 24 for c in centers}))
    print("Center fiber multiplicity gcd:", gcd(*centers.values()))
    print("Centers:", sorted(centers))
    print("Missing odd residues modulo 48:", sorted(set(range(1, 48, 2)) - {v % 48 for v in h}))
    # Exact set cover on the 12 even center phases, with a physical edge as a set.
    cover = {0: ()}
    for j in range(len(EDGES)):
        phases = {((h[s] + d[s][j] // 2) % 24) // 2 for s in range(len(h))
                  if abs(d[s][j]) == 2}
        mask = sum(1 << x for x in phases)
        for seen, old in list(cover.items()):
            new = seen | mask
            if new not in cover or len(old) + 1 < len(cover[new]):
                cover[new] = old + (j,)
    portfolio = cover[(1 << 12) - 1]
    print("Minimum number of marked arcs covering all 12 center phases:", len(portfolio))
    for j in portfolio:
        marks = sorted({(h[s] + d[s][j] // 2) % 24 for s in range(len(h)) if abs(d[s][j]) == 2})
        print("  edge", EDGES[j], "center phases", marks)
    for s, j in ((0, 0), (79, 8), (255, 35)):
        assert values_at(CODE ^ (1 << j), s << 1, False) == h[s] + d[s][j]
    j = EDGES.index((1, 2))
    cycle_seams = [values_at(CODE ^ (1 << j), s << 1, cycles=True)
                   - values_at(CODE, s << 1, cycles=True) for s in range(len(h))]
    coupled = transform([a - 6*b for a, b in zip(tables[j], cycle_seams)], -1)
    edge_mask = 3
    assert all(v == 0 for s, v in enumerate(coupled)
               if s.bit_count() > 4 or
               (s.bit_count() == 4 and (s & edge_mask).bit_count() != 1))
    visible = sum(s.bit_count() == 4 and bool(v) for s, v in enumerate(coupled))
    print("Quartic terms of the marked H-6C potential for edge (1,2):", visible,
          "all incident at exactly one endpoint")
    # Select 182 cuts using finite-field pivots. Nonzero determinant modulo
    # this auxiliary prime certifies invertibility over the rationals.
    e = {1, 2}
    characters = [w for w in range(1, 1 << N)
                  if w.bit_count() in (2, 4, 6) and not e <= set(bits(w))
                  and (w.bit_count() != 6 or e & set(bits(w)))]
    assert len(characters) == 182
    pivot = {}
    p = 1009
    for s in range(1 << (N - 1)):
        row = [(-1) ** ((w & (s << 1)).bit_count()) for w in characters]
        for j, basis in pivot.items():
            f = row[j] % p
            if f:
                row = [(v - f * b) % p for v, b in zip(row, basis)]
        j = next((j for j, v in enumerate(row) if v % p), None)
        if j is not None:
            pivot[j] = [(v * pow(row[j], -1, p)) % p for v in row]
    assert len(pivot) == len(characters)
    print("Exact linear seam observation dimension:", len(pivot))


if __name__ == "__main__":
    main()
