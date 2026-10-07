#!/usr/bin/env python3
"""Audit A: completeness of the session's +1 head grammar (heads_plus1.py: debt 3, |v| = |u|+3, sum u = sum v + 2,
F_v(1) = 4 F_u(1) + 1, i.e. the child reaches 4a+1 when the source is at a).

Own method: enumerate EVERY absorbing post-run parity pattern of the D = 3 chain from the universal state (3, -26)
(depth-first, exact integer chain of a6_numerics.step) for r = 1 (prefix 1,1) and r = 2 (prefix 1,0,0) up to depth LMAX.
For each pattern (first absorption at s) rebuild Y mod 2^s, X = 27 Y - 26, and read the merge shape from the actual
parities: exactly one of the two orbits is odd at time s-1; if it is the source, the child's last odd point b' satisfies
b' = 4^i a + (4^i - 1)/3 (child-side ladder of index i, final letters c, c+2i); otherwise the source is on the ladder.
The session's grammar is the child-side ladder with i = 1 only.  We report the share of each shape and check the
general head identity F_v(1) = 4^i F_u(1) + (4^i - 1)/3 (resp. with u, v swapped) on every pattern.
"""
import sys
from collections import Counter
from fractions import Fraction as Fr
from a6_numerics import step

LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 19

def T(x):
    return (3*x + 1) >> 1 if x & 1 else x >> 1

def dfs(prefix):
    out = []
    stack = [((), 3, -26)]
    while stack:
        bits, k, E = stack.pop()
        if len(bits) >= LMAX: continue
        nb_choices = (prefix[len(bits)],) if len(bits) < len(prefix) else (0, 1)
        for b in nb_choices:
            k2, E2 = step(k, E, b)
            nb = bits + (b,)
            if k2 == 0 and E2 == 0: out.append(nb)
            else: stack.append((nb, k2, E2))
    return out

def Y_from_bits(bits):
    """invert the Terras parity map: Y mod 2^L with the given parity vector (bit by bit lifting)"""
    Y = 0
    for L in range(1, len(bits) + 1):
        for cand in (Y, Y + (1 << (L - 1))):
            z = cand; ok = True
            for b in bits[:L]:
                if (z & 1) != b: ok = False; break
                z = T(z)
            if ok: Y = cand; break
        else:
            raise RuntimeError("no lift")
    return Y

def letters(par):
    w = []; i = 0
    while i < len(par):
        assert par[i] == 1
        j = i + 1
        while j < len(par) and par[j] == 0: j += 1
        w.append(j - i); i = j
    return w

def Fw(w, x):
    x = Fr(x)
    for a in w: x = (3*x + 1) / Fr(2**a)
    return x

def classify(bits):
    s = len(bits)
    Y = Y_from_bits(bits) + (1 << s) * (2**300 + 1)
    X = 27*Y - 26
    px, py = [], []
    u, v = X, Y
    for _ in range(s):
        px.append(u & 1); py.append(v & 1); u, v = T(u), T(v)
    assert u == v
    if px[-1] == 1 and py[-1] == 0:
        tc = max(t for t in range(s) if py[t] == 1)
        i2 = s - 1 - tc
        assert i2 % 2 == 0
        i = i2 // 2
        uw = letters(px[:s - 1]); vw = letters(py[:tc])     # heads (complete letters before the last odd points)
        ok = (len(vw) == len(uw) + 3 and sum(uw) == sum(vw) + 2*i and
              Fw(vw, 1) == 4**i * Fw(uw, 1) + Fr(4**i - 1, 3))
        return ('child-ladder', i, ok, tuple(uw), tuple(vw))
    else:
        ts = max(t for t in range(s) if px[t] == 1)
        i2 = s - 1 - ts
        assert i2 % 2 == 0
        i = i2 // 2
        uw = letters(px[:ts]); vw = letters(py[:s - 1])
        ok = (len(vw) == len(uw) + 3 and sum(vw) == sum(uw) + 2*i and
              Fw(uw, 1) == 4**i * Fw(vw, 1) + Fr(4**i - 1, 3))
        return ('source-ladder', i, ok, tuple(uw), tuple(vw))

if __name__ == "__main__":
    for r, pre in ((1, (1, 1)), (2, (1, 0, 0))):
        pats = dfs(pre)
        shape_mass = Counter(); shape_cnt = Counter(); allok = True; examples = {}
        for p in pats:
            sh, i, ok, uw, vw = classify(p)
            allok &= ok
            m = 2.0**(-(len(p) - len(pre)))
            shape_mass[(sh, i)] += m; shape_cnt[(sh, i)] += 1
            if (sh, i) not in examples or len(p) < examples[(sh, i)][0]:
                examples[(sh, i)] = (len(p), uw, vw)
        tot = sum(shape_mass.values())
        sess = shape_mass[('child-ladder', 1)]
        print(f"r={r}: {len(pats)} absorbing patterns up to post-run depth {LMAX}; total conditional mass {tot:.6f}; "
              f"general head identity holds on all: {allok}")
        for key in sorted(shape_mass):
            print(f"   {key[0]:13s} i={key[1]}: {shape_cnt[key]:5d} patterns, mass {shape_mass[key]:.6f} "
                  f"({shape_mass[key]/tot:.1%}); shallowest example (depth, u, v) = {examples[key]}")
        print(f"   session grammar (child-ladder, i=1) captures {sess/tot:.1%} of the absorbing mass up to depth {LMAX}")

    # generalized compiler as an identity of affine maps of n (K >= 4, J >= 3), for the shallowest example of each shape
    import random
    rnd = random.Random(8)
    def ident(u, v, i, side, K, J, c):
        if side == 'child-ladder':
            src = (1,)*(K-1) + (2,)*J + u + (c,);        chd = (1,)*(K-4) + (4,1,1) + (2,)*(J-3) + v + (c + 2*i,)
        else:
            src = (1,)*(K-1) + (2,)*J + u + (c + 2*i,);  chd = (1,)*(K-4) + (4,1,1) + (2,)*(J-3) + v + (c,)
        return all(Fw(chd, (n + 1)/8 - 1) == Fw(src, n) for n in (Fr(0), Fr(1)))
    cnt = 0
    for r, pre in ((1, (1, 1)), (2, (1, 0, 0))):
        for p in dfs(pre):
            sh, i, ok, uw, vw = classify(p)
            K = rnd.randint(4, 40); J = rnd.randint(3, 30); c = rnd.randint(1, 5)
            assert ident(uw, vw, i, sh, K, J, c), (sh, i, uw, vw)
            cnt += 1
    print(f"generalized compiler (ladder index i, either side) is an identity of affine maps for all {cnt} patterns (random K, J, c)")
