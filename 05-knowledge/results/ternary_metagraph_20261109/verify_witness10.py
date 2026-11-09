"""Independent re-check of the order-10 witnesses written by ternary_order10.c in witness mode (lines 'H_3 code').

For each witness: decode the 45-bit code (pairs i < j in lexicographic order, bit set iff i -> j), check that the
tournament is strong, count the directed cycles on every odd vertex subset (Hamiltonian-cycle DP from the least
vertex, written independently of the C code), form the packing numbers alpha_1, alpha_2, alpha_3 by summing over
disjoint subsets (at order 10 at most three odd cycles are pairwise disjoint), and check
    H_3 = 1 + 3 alpha_1 + 9 alpha_2 + 27 alpha_3   equals the claimed value,
    H   = 1 + 2 alpha_1 + 4 alpha_2 + 8 alpha_3    equals the Hamiltonian-path count (OCF, independent DP).
usage: python verify_witness10.py WIT_GLOB [sample_size]   (sample_size 0 = all)
"""
import sys
import glob
import random
from multiprocessing import Pool

N = 10


def decode(code):
    out = [0] * N
    b = 0
    for i in range(N):
        for j in range(i + 1, N):
            if code >> b & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1
    return out


def strong(out):
    full = (1 << N) - 1
    reach = 1
    while True:
        nxt = reach
        for v in range(N):
            if reach >> v & 1:
                nxt |= out[v]
        if nxt == reach:
            break
        reach = nxt
    if reach != full:
        return False
    inn = [0] * N
    for v in range(N):
        for u in range(N):
            if out[v] >> u & 1:
                inn[u] |= 1 << v
    reach = 1
    while True:
        nxt = reach
        for v in range(N):
            if reach >> v & 1:
                nxt |= inn[v]
        if nxt == reach:
            break
        reach = nxt
    return reach == full


def cycles_by_subset(out):
    """cyc[S] = number of directed Hamiltonian cycles of T[S] for odd |S| >= 3"""
    cyc = {}
    for s in range(N):
        # paths starting at s that use only vertices > s
        f = {(1 << s, s): 1}
        frontier = [(1 << s, s)]
        layer = {(1 << s, s): 1}
        for _ in range(N - s - 1):
            nl = {}
            for (S, v), c in layer.items():
                nb = out[v] & ~S & ~((1 << (s + 1)) - 1)
                while nb:
                    u = (nb & -nb).bit_length() - 1
                    nb &= nb - 1
                    key = (S | (1 << u), u)
                    nl[key] = nl.get(key, 0) + c
            for key, c in nl.items():
                f[key] = f.get(key, 0) + c
            layer = nl
        for (S, v), c in f.items():
            k = bin(S).count("1")
            if k >= 3 and k % 2 == 1 and out[v] >> s & 1:
                cyc[S] = cyc.get(S, 0) + c
    return cyc


def hpaths(out):
    f = [[0] * N for _ in range(1 << N)]
    for v in range(N):
        f[1 << v][v] = 1
    for S in range(1, 1 << N):
        row = f[S]
        for v in range(N):
            c = row[v]
            if c:
                nb = out[v] & ~S
                while nb:
                    u = (nb & -nb).bit_length() - 1
                    nb &= nb - 1
                    f[S | (1 << u)][u] += c
    return sum(f[(1 << N) - 1])


def check(item):
    val, code = item
    out = decode(code)
    if not strong(out):
        return (val, code, "not strong")
    cyc = cycles_by_subset(out)
    keys = list(cyc)
    a1 = sum(cyc.values())
    a2 = 0
    a3 = 0
    tri = [S for S in keys if bin(S).count("1") == 3]
    for i, S in enumerate(keys):
        for T in keys[i + 1:]:
            if S & T == 0:
                a2 += cyc[S] * cyc[T]
    for i, S in enumerate(tri):
        for j in range(i + 1, len(tri)):
            T = tri[j]
            if S & T:
                continue
            for U in tri[j + 1:]:
                if U & (S | T) == 0:
                    a3 += cyc[S] * cyc[T] * cyc[U]
    h3 = 1 + 3 * a1 + 9 * a2 + 27 * a3
    if h3 != val:
        return (val, code, f"H_3 {h3}")
    if 1 + 2 * a1 + 4 * a2 + 8 * a3 != hpaths(out):
        return (val, code, "OCF fails")
    return None


def main():
    items = {}
    for f in sorted(glob.glob(sys.argv[1])):
        for line in open(f):
            v, c = line.split()
            items.setdefault(int(v), int(c))
    vals = sorted(items)
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    if k and k < len(vals):
        vals = sorted(random.Random(87).sample(vals, k))
    with Pool() as pool:
        res = pool.map(check, [(v, items[v]) for v in vals], chunksize=8)
    bad = [r for r in res if r]
    print(f"witnesses: {len(items)} distinct values; checked {len(vals)}; failures {len(bad)} {bad[:5]}")


if __name__ == "__main__":
    main()
