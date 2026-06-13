"""SIGNED LRC homometry -- anatomy of collisions + Galois-stability of the zero set.
monad-explorer-2026-06-06-S705.  Supports THM-415 (prime => full orbit) and HYP-2273 (count law).

For composite C: list every collision pair, the Diff set (flipped runners), and classify by
gcd(i,C) content.  Verify the sine-transform M is invertible and that the zero set
{t: Phi(delta)_t=0} is closed under t->a*t (Galois) -- the engine of THM-415.
"""
from itertools import product
import math
from fractions import Fraction
import numpy as np


def folded_spectrum(V, eps, C):
    out = []
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            f = (eps[i] * V[i] - eps[j] * V[j]) % C
            out.append(min(f, C - f))
    return tuple(sorted(out))


def phi(eps, V, C):
    """signed sine vector Phi(eps)_t = sum_i eps_i sin(2 pi t i / C), t=1..(C-1)//2."""
    return [sum(eps[i] * math.sin(2 * math.pi * t * V[i] / C) for i in range(len(V)))
            for t in range(1, (C - 1) // 2 + 1)]


def sine_matrix(C):
    h = (C - 1) // 2
    return np.array([[math.sin(2 * math.pi * t * i / C) for i in range(1, h + 1)]
                     for t in range(1, h + 1)])


print("M = sine-transform invertibility (det should be nonzero):")
for n in range(4, 16):
    C = 2 * n - 1
    M = sine_matrix(C)
    print(f"  C={C:>3}: |det| = {abs(np.linalg.det(M)):.3e}  (h={(C-1)//2})")

print("\nGALOIS-STABILITY test of zero set {t: Phi(delta)_t=0} under t->a*t (a in (Z/C)*):")
print("  (for random delta supported on a subset Diff; zero set must be a union of gcd-classes)")
for C in [21, 25, 27, 35]:
    h = (C - 1) // 2
    V = list(range(1, h + 1))
    # random-ish delta
    rng = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1]
    Diff = list(range(1, min(h, 6) + 1))
    delta = [0] * h
    for k, i in enumerate(Diff):
        delta[i - 1] = rng[k % len(rng)]
    ph = [sum(delta[i - 1] * math.sin(2 * math.pi * t * i / C) for i in range(1, h + 1))
          for t in range(1, h + 1)]
    zero = set(t for t in range(1, h + 1) if abs(ph[t - 1]) < 1e-9)
    # express zero set in Z/C symmetric form and check closure under multiplication
    zfull = set()
    for t in zero:
        zfull.add(t % C); zfull.add((-t) % C)
    units = [a for a in range(1, C) if math.gcd(a, C) == 1]
    closed = all(((a * t) % C in zfull) for t in zfull for a in units)
    gcds = sorted(set(math.gcd(t, C) for t in zero)) if zero else []
    print(f"  C={C:>3}: zero-set (in 1..{h}) = {sorted(zero)}  gcd-classes={gcds}  Galois-closed={closed}")

print("\nANATOMY: every collision pair, Diff set, and gcd(Diff,C) content:")
for n in range(5, 19):
    C = 2 * n - 1
    if all(C % d for d in range(2, int(C**0.5) + 1)):  # prime
        continue
    if n > 18:
        break
    V = list(range(1, n))
    fold_map = {}
    for bits in product([0, 1], repeat=n - 2):
        eps = tuple([1] + [1 if b else -1 for b in bits])
        fold_map.setdefault(folded_spectrum(V, list(eps), C), []).append(eps)
    coll = [g for g in fold_map.values() if len(g) > 1]
    print(f"\n  C={C} (n={n}): {len(coll)} collision classes, deficiency={sum(len(g)-1 for g in coll)}")
    shown = 0
    diffstats = {}
    for g in coll:
        base = g[0]
        for other in g[1:]:
            Diff = [V[i] for i in range(len(V)) if base[i] != other[i]]
            Diffc = [V[i] for i in range(len(V)) if (base[i] != other[i]) ^ True]
            Diff = min([Diff, [V[i] for i in range(len(V)) if base[i] == other[i]]], key=len)
            gc = tuple(sorted(set(math.gcd(d, C) for d in Diff)))
            key = (len(Diff), gc)
            diffstats[key] = diffstats.get(key, 0) + 1
            if shown < 6:
                print(f"      Diff(flip)={Diff}  |Diff|={len(Diff)}  gcd-content={gc}")
                shown += 1
    print(f"      Diff-stats {{(|Diff|,gcd-set): count}} = {diffstats}")
