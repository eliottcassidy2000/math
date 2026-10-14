#!/usr/bin/env python3
"""Audit A (independent), THM-4640 (3): the uniform switches induced by primitive collisions, on actual integers.

For every canonical primitive pair (s, t), |s| < |t|, D = |t| - |s|, exponents <= E:
  main switch: sources n whose U-word starts 1^r (head(s), c), r >= D;   m = (n+1)/2^D - 1
  twin switch: sources n whose U-word starts 1^r (2, e_1 - 1, ..., c) (root twin (1).s), r >= D - 1;
               m = (n+1)/2^(D-1) - 1   (reset-2 sources: the letter after the ones-run is 2)
We check U^j(n) = U^j(m) with j = r + |shorter word|, that m < n, and that m's word starts 1^(r - D') (longer word).
Integers with a prescribed U-word prefix are built from the Terras parity vector (bijective mod 2^L).
Usage: python3 a6_switch_integers.py [E=14] [NSRC=3]"""
import sys, itertools, random, time
from collections import defaultdict

E = int(sys.argv[1]) if len(sys.argv) > 1 else 14
NSRC = int(sys.argv[2]) if len(sys.argv) > 2 else 3
random.seed(20261014)
t0 = time.time()


def K_of(s):
    K = -1
    for e in s:
        K = 3 * K + (1 << e)
    return K


def head_word(rep):
    """head word with increasing representation rep: letters e_1 + 1, e_2 - e_1, ..."""
    if not rep:
        return ()
    return (rep[0] + 1,) + tuple(rep[j] - rep[j - 1] for j in range(1, len(rep)))


def U(n):
    n = 3 * n + 1
    while n % 2 == 0:
        n //= 2
    return n


def Uj(n, j):
    for _ in range(j):
        n = U(n)
    return n


def word(n, j):
    out = []
    for _ in range(j):
        n = 3 * n + 1
        a = 0
        while n % 2 == 0:
            n //= 2; a += 1
        out.append(a)
    return tuple(out)


def terras_parities(x, L):
    out = []
    for _ in range(L):
        out.append(x & 1)
        x = (3 * x + 1) // 2 if x & 1 else x // 2
    return out


def prefix_class(w):
    """residue of odd n mod 2^(A+1) whose U-word starts with exactly w"""
    pv = []
    for a in w:
        pv += [1] + [0] * (a - 1)
    pv.append(1)
    n = 0
    for k in range(len(pv)):
        for cand in (n, n + (1 << k)):
            if terras_parities(cand + (1 << (len(pv) + 2)), k + 1) == pv[:k + 1]:   # positive representative
                n = cand
                break
        else:
            raise RuntimeError("no lift")
    return n, 1 << len(pv)


reps = defaultdict(list)
for nn in range(E + 1):
    for s in itertools.combinations(range(2, E + 1), nn):
        reps[K_of(s)].append(s)
fam = 0; srcs = 0
for k, rs in reps.items():
    if k == -1:
        continue
    for s, t in itertools.combinations(rs, 2):
        if s[-1] == t[-1]:
            continue
        if len(s) > len(t):
            s, t = t, s
        D = len(t) - len(s)
        cmin = max(1, 1 + t[-1] - s[-1])
        for c in (cmin, cmin + 1, cmin + 3):
            c2 = c + s[-1] - t[-1]
            short, long_ = head_word(s) + (c,), head_word(t) + (c2,)
            twin = head_word((1,) + s) + (c,)
            assert twin[0] == 2
            for (sh, lg, Dsw, kind) in ((short, long_, D, "main"), (twin, long_, D - 1, "twin")):
                assert Dsw >= 1
                for r in (Dsw, Dsw + 1, Dsw + 3):
                    w = (1,) * r + sh
                    res, mod = prefix_class(w)
                    for _ in range(NSRC):
                        n = res + mod * random.randrange(1, 10 ** 6)
                        assert word(n, len(w)) == w
                        mm = (n + 1) // 2 ** Dsw - 1
                        assert (n + 1) % 2 ** Dsw == 0 and 0 < mm < n
                        j = r + len(sh)
                        assert Uj(n, j) == Uj(mm, j), (kind, sh, lg, r, n)
                        assert word(mm, j) == (1,) * (r - Dsw) + lg, (kind, sh, lg, r, n, word(mm, j))
                        srcs += 1
                    fam += 1
print(f"E <= {E}: {fam} (pair, c, main/twin, r) families, {srcs} integer sources: every source merges with "
      f"(n+1)/2^D' - 1 at U-step r + |shorter|, and the partner's word is 1^(r-D') (longer word); "
      f"D' = D (main) or D - 1 (twin)  [{time.time()-t0:.1f}s]")
# THM-4555 (vi) example
for c in range(1, 6):
    assert K_of((1, 7)) == K_of((3, 4, 5)) == 125
    assert head_word((1, 7)) + (c,) == (2, 6, c) and head_word((3, 4, 5)) + (c + 2,) == (4, 1, 1, c + 2)
print("THM-4555 (vi): (2,6,c) ~ (4,1,1,c+2) is the root twin of (7) ~ (3,4,5): OK")
