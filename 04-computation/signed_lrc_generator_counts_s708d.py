"""Exact per-generator silent counts A(D) and co-silence (size-4) structure  (monad-explorer-S708).

Method: hash every cut eps (transversal eps_0=+1) to a class id via Phi^2 (validated exact for
C<=39 in s708/s708c). Cut index `bits` has bit (x-2) = sign of magnitude x (x=2..n-1); magnitude 1
fixed +. A flip pattern D (subset of magnitudes, never containing 1) -> Dmask = OR(1<<(x-2)).
Then A(D) = #{ eps : class(bits) == class(bits ^ Dmask) }, exact given the hashing.
Co-silent pairs (D1,D2) with class(bits)==class(bits^D1)==class(bits^D2) build the size-4 classes.
"""
import numpy as np, math
from collections import Counter


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


def Vspan(C):
    halves = subgroup_halves(C)
    half_n = (C - 1) // 2

    def mask(s):
        m = 0
        for x in s:
            m |= 1 << (x - 1)
        return m
    gens = [(d, mask(h)) for d, h in sorted(halves.items())]
    elems = {0: ()}
    for d, g in gens:
        for e, lab in list(elems.items()):
            ne = e ^ g
            if ne not in elems:
                elems[ne] = lab + (d,)
    return elems, halves


def classids(n, ndec=7, batch=1 << 18):
    C = 2 * n - 1
    n1 = n - 1
    Hh = np.arange(1, n1 + 1)
    S = np.sin(2 * math.pi * np.outer(Hh, Hh) / C)
    nfree = n1 - 1
    total = 1 << nfree
    ids = np.empty(total, dtype=np.int32)
    table = {}
    nxt = 0
    for start in range(0, total, batch):
        end = min(start + batch, total)
        idx = np.arange(start, end, dtype=np.uint64)
        eps = np.ones((n1, end - start), dtype=np.int8)
        for b in range(nfree):
            bit = ((idx >> np.uint64(b)) & np.uint64(1)).astype(np.int8)
            eps[b + 1, :] = np.where(bit == 1, 1, -1)
        Phi2 = np.round((S @ eps) ** 2, ndec)
        for j, row in enumerate(Phi2.T):
            key = row.tobytes()
            cid = table.get(key)
            if cid is None:
                cid = nxt; table[key] = cid; nxt += 1
            ids[start + j] = cid
    return C, ids, nxt


def vmask_to_bits(Vmask):
    """convert span-mask (bit x-1 for magnitude x) to cut-index mask (bit x-2, magnitude>=2)."""
    out = 0
    x = 1
    m = Vmask
    while m:
        if m & 1:
            assert x >= 2, "flip pattern unexpectedly contains magnitude 1"
            out |= 1 << (x - 2)
        m >>= 1
        x += 1
    return out


def analyze(n):
    elems, halves = Vspan(2 * n - 1)
    C, ids, ncl = classids(n)
    total = len(ids)
    Vmasks = [m for m in elems if m]
    print(f"\n=== C={C} = {'x'.join(map(str,factor(C)))}  n={n}  |V|={len(elems)} ===")
    # A(D)
    bitmasks = {}
    A = {}
    allidx = np.arange(total, dtype=np.int64)
    for D in Vmasks:
        bm = vmask_to_bits(D)
        bitmasks[D] = bm
        A[D] = int(np.count_nonzero(ids == ids[allidx ^ bm]))
    for D in sorted(Vmasks, key=lambda m: (bin(m).count('1'), m)):
        mags = [i + 1 for i in range(C) if (D >> i) & 1]
        print(f"   D={str(mags):<26} = H_{elems[D]:}  A(D)={A[D]}")
    # size-4 structure: for each eps, collect silent D's; record the GROUP (frozenset) per eps
    silent_sets = Counter()
    # vectorized: for each eps, which D silent? build boolean matrix lazily (|V| small)
    silentbool = {D: (ids == ids[allidx ^ bitmasks[D]]) for D in Vmasks}
    # determine G_eps as frozenset of silent D-masks (plus 0)
    sizes = np.ones(total, dtype=np.int64)
    for D in Vmasks:
        sizes += silentbool[D].astype(np.int64)   # |G_eps| = 1 + #silent generators (only if group!)
    sizehist = Counter(sizes.tolist())
    # report size-4 generator combos
    big = np.where(sizes >= 4)[0]
    combo = Counter()
    for i in big[:200000]:
        G = tuple(sorted(elems[D] for D in Vmasks if silentbool[D][i]))
        combo[G] += 1
    classhist = {s: c // s for s, c in sorted(sizehist.items())}
    defic = sum((s - 1) * c for s, c in classhist.items())
    print(f"   |G_eps| histogram (size:#eps) = {dict(sorted(sizehist.items()))}")
    print(f"   class hist (size:#classes)    = {classhist}   deficiency={defic}")
    if combo:
        print(f"   size>=4 generator-combos (decomp tuples -> #eps): {dict(combo)}")
    return C, A, classhist, defic, elems


if __name__ == "__main__":
    for n in [14, 23]:    # C=27 (3^3), C=45 (3^2 x 5)
        analyze(n)
