#!/usr/bin/env python3
"""Wave 23, P47: search for Jacobi-symbol invariants (Rickards-Stange type) on tournament transfer orbits.

Orbits: states (H,E,L,h) of O/Z/J root-change words from (1,1,1,1) (the GL3-dense, non-abelian system), and,
for contrast, of one-apex words (letters O,Z only; triangular).  For every pair of primitive integer linear
forms a.v, b.v (coefficients in [-1,1], 4 coordinates), over orbit states with b.v odd and positive and
gcd(a.v,b.v)=1, record the frequency of Jacobi symbol +1.  An invariant (frequency 0 or 1) on a non-congruence
pair would be a reciprocity obstruction in the sense of Rickards-Stange / Haag-Kertzer-Rickards-Stange.
"""
import itertools
import random
import sys
from math import gcd


def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    res = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            res = -res
        a %= n
    return res if n == 1 else 0


def step(s, m):
    H, E, L, h = s
    if m == 'O':
        return (H, 0, 2 * L, h)
    if m == 'Z':
        return (H + L - E, h, 2 * L - h, h)
    return (H + L - E, H - E, H + L, H)


def orbit(letters, depth, cap):
    seen = set()
    frontier = [(1, 1, 1, 1)]
    seen.add(frontier[0])
    for _ in range(depth):
        nxt = []
        for s in frontier:
            for m in letters:
                t = step(s, m)
                if t not in seen:
                    seen.add(t)
                    nxt.append(t)
        frontier = nxt
        if len(seen) > cap:
            break
    return list(seen)


def forms():
    out = []
    for c in itertools.product((-1, 0, 1), repeat=4):
        if any(c):
            # primitive up to sign: first nonzero coefficient positive
            first = next(x for x in c if x)
            if first > 0:
                out.append(c)
    return out


def analyse(name, states, nsample=40000, minuse=2000, nondegenerate=False):
    random.seed(7)
    S = random.sample(states, min(nsample, len(states)))
    F = forms()
    vals = {f: [sum(a * b for a, b in zip(f, s)) for s in S] for f in F}
    flagged = []
    npairs = 0
    for fa in F:
        A = vals[fa]
        for fb in F:
            if fa == fb:
                continue
            B = vals[fb]
            plus = use = 0
            for a, b in zip(A, B):
                ok = b > 0 and b % 2 == 1 and gcd(a, b) == 1
                if nondegenerate:
                    ok = ok and b >= 3 and abs(a) >= 3
                if ok:
                    use += 1
                    if jacobi(a, b) == 1:
                        plus += 1
            if use >= minuse:
                npairs += 1
                fr = plus / use
                analyse.lo = min(getattr(analyse, 'lo', 1.0), fr)
                analyse.hi = max(getattr(analyse, 'hi', 0.0), fr)
                if not 0.40 <= fr <= 0.60:
                    flagged.append((fr, use, fa, fb))
    tag = " (values |a.v| <= 2 and b.v = 1 excluded)" if nondegenerate else ""
    print(f"== {name}{tag}: {len(states)} distinct states, sample {len(S)}, pairs with >= {minuse} usable samples: "
          f"{npairs}, outside [0.40,0.60]: {len(flagged)}")
    for fr, use, fa, fb in sorted(flagged)[:60]:
        print(f"   ({fa}.v / {fb}.v): freq(+1) = {fr:.4f} over {use}")
    print(f"   frequency range over all {npairs} pairs: [{analyse.lo:.4f}, {analyse.hi:.4f}]")
    analyse.lo, analyse.hi = 1.0, 0.0
    return npairs, flagged


if __name__ == '__main__':
    ozj = orbit('OZJ', 14, 3_000_000)
    oz = orbit('OZ', 22, 3_000_000)
    analyse("O/Z/J root-change orbit (non-abelian, blocks GL3-dense)", ozj)
    analyse("O/Z/J root-change orbit (non-abelian, blocks GL3-dense)", ozj, nondegenerate=True)
    analyse("one-apex O/Z orbit (triangular)", oz)
