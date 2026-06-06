"""Silent-move RANK law for signed-LRC AP_n deficiency.  monad-explorer-S708.

Generators: G = {H_d : d|C, 1<d<C}, H_d = half-system of subgroup K_d=<C/d> (speeds in 1..n1).
For each cut eps, SILENT(eps) = {g in <G> : flipping g preserves Phi^2}.  Claim: SILENT(eps) is a
subGROUP; class-size(eps)=|SILENT(eps)|=2^{rank(eps)}.  deficiency = sum_eps (1 - 2^{-rank(eps)}).

We (1) verify class-size == 2^{rank} where rank = dim_GF2 span of silent generators (realized & tested),
(2) tabulate the rank distribution (#cuts at each rank) -> the clean count-law ingredients,
(3) print N(g) for every nontrivial group element g (= #cuts where g silent).
"""
import numpy as np
import math
import sys
from collections import Counter, defaultdict
from itertools import combinations


def factor_str(C):
    f, m, d = [], C, 2
    while d * d <= m:
        while m % d == 0:
            f.append(d); m //= d
        d += 1
    if m > 1:
        f.append(m)
    return "x".join(map(str, f))


def proper_divisors(C):
    return [d for d in range(2, C) if C % d == 0]


def halfmask(C, d, n1):
    g = C // d
    Kd = set((g * j) % C for j in range(d))
    mask = np.zeros(n1, dtype=bool)
    for x in Kd:
        if 0 < x <= n1:
            mask[x - 1] = True
    return mask


def analyze(C):
    n1 = (C - 1) // 2
    H = np.arange(1, n1 + 1)
    S = np.sin(2 * math.pi * np.outer(H, H) / C)
    nfree = n1 - 1
    total = 1 << nfree
    divs = proper_divisors(C)
    gmasks = {d: halfmask(C, d, n1) for d in divs}
    # build cuts
    eps = np.ones((n1, total), dtype=np.int8)
    idx = np.arange(total, dtype=np.uint64)
    for b in range(nfree):
        bit = ((idx >> np.uint64(b)) & np.uint64(1)).astype(np.int8)
        eps[b + 1, :] = np.where(bit == 1, 1, -1)
    Phi2 = np.round((S @ eps) ** 2, 6)
    # class sizes (ground truth)
    cls = defaultdict(int)
    keys = [Phi2[:, j].tobytes() for j in range(total)]
    for k in keys:
        cls[k] += 1
    class_size = np.array([cls[keys[j]] for j in range(total)])
    n_classes = len(cls)
    defic = total - n_classes

    # FULL GROUP <H_d>: enumerate XOR of every subset of generator masks, canonicalize mod global
    # complement, keep nontrivial distinct elements. Test EACH for silence per cut (a combined move
    # can be silent while neither factor is -- so we must test the whole group, not just generators).
    full = np.ones(n1, dtype=bool)
    elems = {}
    ndiv = len(divs)
    for sub in range(1, 1 << ndiv):
        v = np.zeros(n1, dtype=bool)
        for r in range(ndiv):
            if sub & (1 << r):
                v ^= gmasks[divs[r]]
        if not v.any():
            continue
        comp = full & ~v
        key = min(v.tobytes(), comp.tobytes())  # canonicalize mod global complement
        elems[key] = v
    group_elems = list(elems.values())   # nontrivial group elements (mod global flip)

    silent_elem = np.zeros((len(group_elems), total), dtype=bool)
    for r, v in enumerate(group_elems):
        flip = np.where(v[:, None], -1, 1).astype(np.int8)
        Phi2f = np.round((S @ (eps * flip)) ** 2, 6)
        same = np.all(np.abs(Phi2f - Phi2) < 1e-4, axis=0)
        changed = np.any((eps * flip) != eps, axis=0)
        silent_elem[r] = same & changed

    # silent subgroup size = (#silent group elements among these) + 1 (identity); rank = log2 size
    nsil = silent_elem.sum(axis=0)             # #nontrivial silent elements per cut
    size_pred = nsil + 1
    rank = np.zeros(total, dtype=np.int64)
    # rank = log2(size) where size is a power of 2; compute safely
    for j in range(total):
        s = int(size_pred[j])
        rank[j] = s.bit_length() - 1 if (s & (s - 1)) == 0 else -1  # -1 flags non-power-of-2

    # verify class_size == 2^rank
    ok = bool(np.all(class_size == (1 << rank)))
    rank_dist = Counter(rank.tolist())
    defic_formula = sum(c * (1 - 2.0 ** (-r)) for r, c in rank_dist.items())

    print(f"=== C={C} ({factor_str(C)}) total={total} defic={defic} "
          f"classes={n_classes}  class_size==2^rank? {ok}")
    print(f"    rank distribution (rank:#cuts) = {dict(sorted(rank_dist.items()))}")
    print(f"    deficiency from rank formula = {defic_formula:.4f}  (matches {abs(defic_formula-defic)<1e-6})")
    # size-k class counts
    sizehist = Counter(cls.values())
    print(f"    class-size hist = {dict(sorted(sizehist.items()))}")
    return defic, dict(sorted(rank_dist.items())), dict(sorted(sizehist.items()))


targets = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [9, 15, 21, 25, 27, 33, 35, 39, 45]
results = {}
for C in targets:
    results[C] = analyze(C)
    print(flush=True)
print("SUMMARY (C: defic, rank_dist, size_hist):")
for C, (d, rd, sh) in results.items():
    print(f"  {C} ({factor_str(C)}): defic={d} rank={rd} sizes={sh}")
