#!/usr/bin/env python3
"""THM-480 (claudebox-2026-06-11-S2): the skew tower's tournament-gauge row code is the
d+ ladder: e8 = RM(1,3) at order 8, d16+ at 16, d32+-enumerator at 32; Sylvester contrast
RM(1,m). Pure python, exact."""
import sys, itertools
from collections import Counter
sys.path.insert(0, '/home/claude/math-research/04-computation')
from tournament_barba_flag_cbx1 import paley_S, skewH_from_drt, double_skewH

def binrows(H): return [int(''.join('1' if x < 0 else '0' for x in row), 2) for row in H]
def basis_of(vecs):
    basis = []
    for v in vecs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v); basis.sort(reverse=True)
    return basis
def all_words(basis):
    words = [0]
    for v in basis: words += [w ^ v for w in words]
    return words
def wd(basis): return Counter(bin(w).count('1') for w in all_words(basis))

def sylvester(m):
    H = [[1]]
    for _ in range(m):
        H = [r + r for r in H] + [r + [-x for x in r] for r in H]
    return H

def dplus(m):
    """d_{2m}+ : even number of 11-pairs + glue 1010...10."""
    pair = lambda i: (0b11 << (2 * i))
    gens = [pair(i) ^ pair(i + 1) for i in range(m - 1)]
    g = 0
    for i in range(m): g |= (1 << (2 * i))
    return basis_of(gens + [g])

def rm(r, m):
    pts = list(itertools.product([0, 1], repeat=m))
    gens = []
    for d in range(r + 1):
        for mo in itertools.combinations(range(m), d):
            v = 0
            for k, p in enumerate(pts):
                val = 1
                for i in mo: val &= p[i]
                if val: v |= (1 << k)
            gens.append(v)
    return basis_of(gens)

def indecomposable(basis, n):
    w4 = [w for w in all_words(basis) if bin(w).count('1') == 4]
    par = list(range(n))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for w in w4:
        cs = [i for i in range(n) if (w >> i) & 1]
        for c in cs[1:]: par[find(cs[0])] = find(c)
    return len({find(i) for i in range(n)}) == 1, len(w4)

def self_orthogonal(basis):
    return all(bin(a & b).count('1') % 2 == 0 for a in basis for b in basis)

if __name__ == '__main__':
    H8 = skewH_from_drt(paley_S(7)); H16 = double_skewH(H8); H32 = double_skewH(H16)
    print("== the tower ladder (tournament gauge M = I + S) ==")
    b8 = basis_of(binrows(H8))
    print(f"order 8 : dim {len(b8)}, wd {dict(sorted(wd(b8).items()))}")
    assert len(b8) == 4 and dict(wd(b8)) == {0: 1, 4: 14, 8: 1}
    # explicit coordinate-permutation witness C(H8) ~ RM(1,3) (Type II [8,4] is unique up to perm)
    import itertools as it
    target = set(all_words(rm(1, 3)))
    w8 = all_words(b8)
    witness = None
    for perm in it.permutations(range(8)):
        ok = True
        for w in w8:
            x = 0
            for i in range(8):
                if (w >> i) & 1: x |= (1 << perm[i])
            if x not in target: ok = False; break
        if ok: witness = perm; break
    assert witness is not None
    print(f"  C(H8) ~ RM(1,3) = e8 via coordinate permutation {witness}  (the r=1 self-dual RM point)")
    b16 = basis_of(binrows(H16))
    ind16, n4 = indecomposable(b16, 16)
    print(f"order 16: dim {len(b16)}, wd {dict(sorted(wd(b16).items()))}, "
          f"self-orth {self_orthogonal(b16)}, indecomposable {ind16}")
    assert len(b16) == 8 and self_orthogonal(b16) and ind16
    assert wd(b16) == wd(dplus(8))
    print("  Type II [16,8], indecomposable => d16+ (NOT e8+e8); matches kps1 pin tower-16 ~ had.16.3")
    b32 = basis_of(binrows(H32))
    print(f"order 32: dim {len(b32)}, wd {dict(sorted(wd(b32).items()))}")
    assert len(b32) == 16 and self_orthogonal(b32)
    assert wd(b32) == wd(dplus(16)), "order-32 enumerator != d32+"
    assert all(k % 4 == 0 for k in wd(b32))
    print("  Type II [32,16] with the d32+ weight enumerator (A4 = 120 = C(16,2)); != RM(2,5) (A4=0)")
    print("\n== Sylvester contrast ==")
    for m in (3, 4, 5):
        bs = basis_of(binrows(sylvester(m)))
        bsa = basis_of(binrows(sylvester(m)) + [(1 << (2 ** m)) - 1])
        print(f"order {2**m}: Sylvester row-span dim {len(bs)} (+all-ones: {len(bsa)}); RM(1,{m}) dim {m+1}")
        assert len(bsa) == m + 1
        assert sorted(all_words(bsa)) == sorted(all_words(rm(1, m)))
    print("  Sylvester + 1 == RM(1,m) AS SETS at m=3,4,5 (the biorthogonal end)")
    print("\nall checks passed")
