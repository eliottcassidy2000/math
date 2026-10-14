#!/usr/bin/env python3
"""Haar coverage of downward uniform switches from primitive collisions at -1 (exact rationals).

Source model: odd 2-adic n; its U-word letters are i.i.d. with P(a) = 2^-a (a >= 1). A word is 1^r u with u reduced.
A primitive collision at -1 (two canonical increasing representations s, t of the same K with different last
exponents, |s| < |t|) yields the word collision (head(s), c) ~ (head(t), c') for every c >= c_min (c - c' = e_t - e_s),
and by THM-4555 (iv) a downward switch for every source whose word starts 1^r (head(s), c), c >= c_min, r >= D,
D = |t| - |s| (even, >= 2): n merges at equal time with (n+1)/2^D - 1 < n.

Coverage(E) = P(the source's word starts 1^r u with u in a cylinder (head(s), >= c_min) of some primitive pair with
max exponent <= E and D <= r).  Computed exactly by a trie over heads, separately for each needed run length.

Audit A (2026-10-14) split: c(E) counts canonical pairs and their root twins but NOT the root collision (a) ~ (2, a-2)
itself, whose switch (THM-4555's reset switch, m = (n-1)/2) covers exactly 1/4 (sources 1^r c, r >= 1, c >= 3).
Canonical cylinders lie inside that domain (first post-run letter >= 3, r >= D >= 2); twin cylinders lie on reset-2
sources (first post-run letter 2).  The script prints the canonical part, the twin part, c(E) = canonical + twin, and
the coverage with the root included, and asserts exactly: twin = 2 * canonical, with-root = 1/4 + twin, c <= 3/8.
Usage: python3 switch_coverage.py E"""
import sys, itertools, time
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, '.')
from repair_algebra import K_of, word_from

def primitive_pairs(E):
    reps = defaultdict(list)
    for n in range(E + 1):
        for s in itertools.combinations(range(2, E + 1), n):      # canonical: no leading exponent 1
            reps[K_of(s)].append(s)
    out = []
    for k, rs in reps.items():
        if k == -1: continue
        for s, t in itertools.combinations(rs, 2):
            es, et = s[-1], t[-1]
            if es == et: continue
            if len(s) == len(t): continue
            if len(s) > len(t): s, t = t, s; es, et = et, es
            D = len(t) - len(s)
            c_min = max(1, 1 + et - es)            # shorter's last letter c, longer's c' = c + es - et >= 1
            out.append((word_from(s, 1)[:-1], c_min, D, s, t))
            # root twin of the shorter side: (1,) + s is another representation of K (word starting with 2, a reset-2
            # source); paired with t it switches downward by D - 1 >= 1 (THM-4555 (vi): (2,6,c) ~ (4,1,1,c+2), D = 1)
            s1 = (1,) + s
            out.append((word_from(s1, 1)[:-1], c_min, D - 1, s1, t))
    return out

def coverage(pairs, rmax=None):
    """P(word = 1^r u, u has a cylinder of a pair with D <= r), exactly: for r >= R0 = max D the cylinder set is
    constant, so the tail sum_{r >= R0} 2^-(r+1) mu = 2^-R0 mu is added in closed form (rmax is ignored)."""
    if not pairs:
        return Fr(0)
    R0 = max(D for _, _, D, _, _ in pairs)
    total = Fr(0)
    for r in range(0, R0 + 1):
        cyl = {}
        for head, cmin, D, s, t in pairs:
            if D <= r:
                cyl[head] = min(cyl.get(head, 10 ** 9), cmin)
        mu = union_measure(cyl)
        if r < R0:
            total += Fr(1, 2 ** (r + 1)) * mu      # P(exactly r ones then a reduced continuation) = 2^-(r+1)
        else:
            total += Fr(1, 2 ** R0) * mu           # all r >= R0 at once
    return total

def union_measure(cyl):
    """measure (under the reduced-word law: first letter a >= 2 w.p. 2^(1-a), then i.i.d. 2^-a) of the union of
    cylinders {words with prefix h, next letter >= c}; computed on the trie of heads"""
    children = defaultdict(set)
    for h in cyl:
        for i in range(len(h)): children[h[:i]].add(h[i])
    memo = {}
    def mu(p):
        if p in memo: return memo[p]
        first = (len(p) == 0)
        thr = cyl.get(p)
        val = Fr(0)
        # letters below the threshold (or all letters if no cylinder at p) that lead to deeper cylinders
        for a in children.get(p, ()):
            if thr is not None and a >= thr: continue
            w = Fr(2, 2 ** a) if first else Fr(1, 2 ** a)
            val += w * mu(p + (a,))
        if thr is not None:
            # P(next letter >= thr): reduced first letter (a >= 2): sum_{a >= max(thr,2)} 2^(1-a) ; else 2^(1-thr)
            if first: val += Fr(2, 2 ** max(thr, 2)) * 2
            else: val += Fr(2, 2 ** thr)
        memo[p] = val
        return val
    return mu(())

if __name__ == '__main__':
    E = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    t0 = time.time()
    ROOT = ((), 3, 1, None, None)                 # reset switch: sources 1^r c, r >= 1, c >= 3, partner (n-1)/2
    for e in range(8, E + 1, 2):
        pairs = primitive_pairs(e)
        canon = [p for p in pairs if p[0][0] >= 3]
        twins = [p for p in pairs if p[0][0] == 2]
        assert len(canon) + len(twins) == len(pairs)
        Ds = sorted({D for _, _, D, _, _ in pairs})
        cov = coverage(pairs, 40)
        cc, ct = coverage(canon, 40), coverage(twins, 40)
        withroot = coverage(pairs + [ROOT], 40)
        assert cov == cc + ct and ct == 2 * cc and withroot == Fr(1, 4) + ct and cov <= Fr(3, 8)
        print(f"max exponent <= {e:2d}: {len(pairs):6d} primitive pairs (switch distances {Ds}); "
              f"c(E) = {float(cov):.6f} = canonical {float(cc):.6f} + twin {float(ct):.6f}; "
              f"twin share of reset-2 mass {float(4 * ct):.4f}; with the root (reset switch) {float(withroot):.6f}"
              f"  [{time.time() - t0:.0f}s]", flush=True)
