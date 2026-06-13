"""The A_t*B_t=0 silent-condition LEMMA + divisor-chain criterion for combined moves.
monad-explorer-S708 (building on the shared S708 deficiency thread; credit concurrent monad-explorer).

LEMMA (to verify): for a move D (subset of runners 1..(C-1)/2), flipping D is silent at cut eps
  <=>  A_t * B_t = 0  for ALL t in 1..(C-1)/2,
where  B_t = sum_{i in D} eps_i sin(2pi t i/C)   (the D-part),
       A_t = sum_{i not in D} eps_i sin(2pi t i/C) (the non-D part).
Proof: Phi(eps^D)_t = A_t - B_t, Phi(eps)_t = A_t + B_t, so Phi'^2=Phi^2 <=> -4 A_t B_t = 0.

We (1) verify the lemma against direct Phi^2 testing for all composite C<=45 and every group move D,
   (2) tabulate which combined moves H_B can EVER be silent (A(D)>0) and test the conjecture:
       H_{d1}+...+{dk} is silent-able only when the multiset {d_i} of subgroup-orders is a chain
       under divisibility (d_1 | d_2 | ... ) OR the FULL generator set.  (refine empirically)
"""
import numpy as np
import math
import sys
from itertools import combinations
from collections import Counter


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


def is_chain(ds):
    s = sorted(ds)
    return all(s[i + 1] % s[i] == 0 for i in range(len(s) - 1))


def analyze(C, verify_lemma=True):
    n1 = (C - 1) // 2
    Hf = np.arange(1, n1 + 1)
    S = np.sin(2 * math.pi * np.outer(Hf, Hf) / C)  # S[t-1,i-1]
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
    Phi = S @ eps
    Phi2 = np.round(Phi ** 2, 6)

    # enumerate group elements with their generator-subset label
    ndiv = len(divs)
    full = np.ones(n1, dtype=bool)
    elem_label = {}
    for sub in range(1, 1 << ndiv):
        v = np.zeros(n1, dtype=bool)
        lab = []
        for r in range(ndiv):
            if sub & (1 << r):
                v ^= gmasks[divs[r]]; lab.append(divs[r])
        if not v.any():
            continue
        comp = full & ~v
        key = min(v.tobytes(), comp.tobytes())
        # keep the label with fewest generators (canonical-ish) — but multiple subsets can give same v;
        # we just need ONE label set; store all labels seen for this key
        elem_label.setdefault(key, (v.copy(), []))
        elem_label[key][1].append(tuple(sorted(lab)))

    lemma_ok = True
    A_of = {}
    for key, (v, labs) in elem_label.items():
        Dmask = v
        flip = np.where(Dmask[:, None], -1, 1).astype(np.int8)
        # ground-truth silence
        Phi2f = np.round((S @ (eps * flip)) ** 2, 6)
        changed = np.any((eps * flip) != eps, axis=0)
        gt = np.all(np.abs(Phi2f - Phi2) < 1e-4, axis=0) & changed
        if verify_lemma:
            # lemma: A_t*B_t=0 all t.  B = S restricted to D columns @ eps_D ; A = Phi - B
            Bd = S[:, Dmask] @ eps[Dmask, :]      # (n1,total)
            Ad = Phi - Bd
            lem = np.all(np.abs(Ad * Bd) < 1e-4, axis=0) & changed
            if not np.array_equal(lem, gt):
                lemma_ok = False
        A_of[key] = (int(gt.sum()), labs, Dmask)

    # report which moves are silent-able and whether their (shortest) label is a divisibility chain
    print(f"=== C={C} ({factor_str(C)}) divisors {divs}  LEMMA A_t*B_t=0 matches Phi^2: {lemma_ok}")
    nonzero, zero = [], []
    for key, (a, labs, Dmask) in A_of.items():
        shortest = min(labs, key=len)
        tag = "CHAIN" if is_chain(shortest) else ("FULLGEN" if len(shortest) == ndiv else "non-chain")
        rec = (shortest, a, tag)
        (nonzero if a > 0 else zero).append(rec)
    print(f"   silent-able moves (A>0):")
    for lab, a, tag in sorted(nonzero, key=lambda r: (-r[1])):
        print(f"      H_{str(lab):<18} A={a:<8} [{tag}]")
    if zero:
        zlabs = sorted([(lab, tag) for lab, a, tag in zero])
        print(f"   NEVER-silent moves (A=0): {[(l, t) for l, t in zlabs]}")
    return lemma_ok, nonzero, zero


targets = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [9, 15, 21, 25, 27, 33, 35, 39, 45]
all_ok = True
for C in targets:
    ok, nz, z = analyze(C)
    all_ok = all_ok and ok
    print(flush=True)
print(f"LEMMA verified for all tested C: {all_ok}")
