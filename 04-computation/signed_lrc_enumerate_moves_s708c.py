"""Enumerate ALL silent moves at composite C (signed-LRC AP_n).  monad-explorer-S708.

For each collision class, list its members; the nontrivial differences (flip-patterns) are the
silent moves. Classify each silent move by the SET of runners it flips (up to global complement),
and by that set's structure mod the subgroups K_d.  Goal: identify what generates the 66 size-2
classes at C=27 and the 8620 at C=45 (beyond the canonical subgroup half-systems).
"""
import numpy as np
import math
from collections import Counter, defaultdict


def factor_str(C):
    f, m, d = [], C, 2
    while d * d <= m:
        while m % d == 0:
            f.append(d); m //= d
        d += 1
    if m > 1:
        f.append(m)
    return "x".join(map(str, f))


def build(C):
    n1 = (C - 1) // 2
    H = np.arange(1, n1 + 1)
    S = np.sin(2 * math.pi * np.outer(H, H) / C)
    nfree = n1 - 1
    total = 1 << nfree
    eps = np.ones((n1, total), dtype=np.int8)
    idx = np.arange(total, dtype=np.uint64)
    for b in range(nfree):
        bit = ((idx >> np.uint64(b)) & np.uint64(1)).astype(np.int8)
        eps[b + 1, :] = np.where(bit == 1, 1, -1)
    Phi2 = np.round((S @ eps) ** 2, 6)
    classes = defaultdict(list)
    for j in range(total):
        classes[Phi2[:, j].tobytes()].append(j)
    return n1, eps, classes, total


def flipset_canonical(eps_a, eps_b, n1):
    """Runners (speeds 1..n1) where eps_a, eps_b differ; canonicalize mod global complement."""
    diff = np.nonzero(eps_a != eps_b)[0]  # 0-based runner indices
    fs = frozenset(int(i + 1) for i in diff)  # speeds
    comp = frozenset(range(1, n1 + 1)) - fs
    return min(fs, comp, key=lambda s: (len(s), sorted(s)))


def analyze(C, show=12):
    n1, eps, classes, total = build(C)
    # collect all silent moves (canonical flip-sets) across all classes of size>1
    move_counter = Counter()
    size_hist = Counter()
    for key, members in classes.items():
        sz = len(members)
        size_hist[sz] += 1
        if sz == 1:
            continue
        base = members[0]
        for m in members[1:]:
            fs = flipset_canonical(eps[:, base], eps[:, m], n1)
            move_counter[fs] += 1
    print(f"=== C={C} ({factor_str(C)}) n1={n1} total={total} class-size hist={dict(sorted(size_hist.items()))}")
    # classify silent moves by |flipset| and by membership in subgroups
    g3 = C // 3 if C % 3 == 0 else None
    # describe each distinct move shape: (size, set of residues, how it sits mod small primes)
    by_size = Counter()
    examples = defaultdict(list)
    for fs, cnt in move_counter.items():
        by_size[len(fs)] += 1
        examples[len(fs)].append((sorted(fs), cnt))
    print(f"    #distinct silent moves = {len(move_counter)}; by flip-set size: {dict(sorted(by_size.items()))}")
    for sz in sorted(examples):
        ex = examples[sz][:show]
        print(f"    |flipset|={sz}: {len(examples[sz])} distinct moves. examples (set, #occurrences):")
        for s, c in ex:
            print(f"        {s}   x{c}")
    return move_counter, n1, C


for C in [9, 15, 21, 25, 27]:
    analyze(C)
    print(flush=True)
