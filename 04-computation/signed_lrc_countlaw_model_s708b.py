"""Count-law model for signed-LRC deficiency via the subgroup lattice  (monad-explorer-S708).

For composite C, the homometry class of a cut eps is a coset eps (+) G_eps, where
  V = GF(2)-span{ H_d : d|C, 1<d<C }      (H_d = positive half of the order-d subgroup)
  G_eps = { D in V : flipping D is silent at eps }  (a subgroup of V; class size = |G_eps|).
deficiency = sum_classes (|class|-1) = sum_eps (1 - 1/|G_eps|).

This script computes, EXACTLY (integer difference-multiset, no floats):
  - V and its dimension,
  - A(D) = #{eps in transversal : D silent}  for each D in V\{0},
  - the joint distribution of G_eps  (which subgroups of V occur, with multiplicity),
  - hence the class-size histogram and deficiency,
to discover the recursion for prime powers (predict C=81) and mixed C.
"""
from itertools import product
from collections import defaultdict, Counter


def factor(C):
    f, m, d = [], C, 2
    while d * d <= m:
        while m % d == 0:
            f.append(d); m //= d
        d += 1
    if m > 1:
        f.append(m)
    return f


def subgroup_halves(C):
    out = {}
    for d in range(2, C):
        if C % d:
            continue
        K = set((C // d * j) % C for j in range(d))
        out[d] = frozenset(x for x in K if 1 <= x <= (C - 1) // 2)
    return out


def span(halves, half_n):
    """GF(2) span of the H_d (as bitmasks over magnitudes 1..half_n). Return set of bitmasks + basis labels."""
    def mask(s):
        m = 0
        for x in s:
            m |= 1 << (x - 1)
        return m
    gens = [(d, mask(h)) for d, h in sorted(halves.items())]
    elems = {0: ()}                       # bitmask -> tuple of d's (one decomposition)
    for d, g in gens:
        new = {}
        for e, lab in elems.items():
            ne = e ^ g
            if ne not in elems and ne not in new:
                new[ne] = lab + (d,)
        elems.update(new)
    return elems, gens


def diff_sig(S, C):
    cnt = [0] * C
    for a in S:
        for b in S:
            if a != b:
                cnt[(a - b) % C] += 1
    return tuple(cnt)


def is_silent(eps, D_mask, runners, C):
    """Flip the magnitudes in D_mask; return True iff homometric to eps."""
    S0 = sorted((eps[i] * runners[i]) % C for i in range(len(runners)))
    eps2 = list(eps)
    for i in range(len(runners)):
        if (D_mask >> (runners[i] - 1)) & 1:
            eps2[i] = -eps2[i]
    S1 = sorted((eps2[i] * runners[i]) % C for i in range(len(runners)))
    return diff_sig(S0, C) == diff_sig(S1, C)


def analyze(n):
    C = 2 * n - 1
    runners = list(range(1, n))
    half_n = n - 1
    halves = subgroup_halves(C)
    elems, gens = span(halves, half_n)
    Vmasks = [m for m in elems if m != 0]
    dimV = len(bin(len(elems) - 1)) if False else (len(elems)).bit_length() - 1
    print(f"\n=== C={C} = {'x'.join(map(str,factor(C)))}  n={n} ===")
    print(f"  subgroups d|C (1<d<C): {sorted(halves)}   |V|={len(elems)} (dim {dimV})")

    # exact enumeration over transversal eps_0=+1
    nfree = n - 2
    A = Counter()                         # D_mask -> #eps silent
    Gdist = Counter()                     # frozenset(silent D masks incl 0) -> count of eps
    for bits in range(1 << nfree):
        eps = [1] * half_n
        for b in range(nfree):
            if not (bits >> b) & 1:
                eps[b + 1] = -1
        # find silent generators among Vmasks
        G = [0]
        for D in Vmasks:
            if is_silent(eps, D, runners, C):
                A[D] += 1
                G.append(D)
        Gdist[frozenset(G)] += 1
    # class size histogram: each eps in a class of size |G_eps|
    sizehist = Counter()
    for G, cnt in Gdist.items():
        sizehist[len(G)] += cnt          # cnt eps's each have class-size len(G)
    # number of classes of each size = (#eps with that size)/size
    classhist = {s: c // s for s, c in sizehist.items()}
    defic = sum((s - 1) * (c) for s, c in classhist.items())
    print(f"  A(D) per generator (D as sorted magnitudes -> count):")
    for D in sorted(Vmasks, key=lambda m: (bin(m).count('1'), m)):
        mags = [i + 1 for i in range(half_n) if (D >> i) & 1]
        lab = elems[D]
        print(f"      {str(mags):<28} decomp H_{lab}   A={A[D]}")
    print(f"  class-size histogram (size -> #classes): {dict(sorted(classhist.items()))}")
    print(f"  deficiency (model, exact) = {defic}")
    return C, classhist, defic, dict(A)


if __name__ == "__main__":
    for n in [5, 14]:        # C=9 (3^2), C=27 (3^3): the prime-power chain
        analyze(n)
    # mixed small: C=15,21,25,35 quick; C=45 too big for this exact per-eps loop (do via histogram script)
    for n in [8, 11, 13, 18]:
        analyze(n)
