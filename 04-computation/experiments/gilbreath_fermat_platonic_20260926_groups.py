#!/usr/bin/env python3
"""gilbreath_fermat_platonic_20260926_groups.py -- the five Platonic solids and the fields of size 3, 4, 5
(session gilbreath-fermat-platonic-20260926, opus, 2026-09-26).

 * Schlafli count: {p, q} with (p-2)(q-2) < 4 (five solids) and = 4 (three plane tilings) [classical; HYP-3772].
 * The rotation groups of the solids are the projective groups over the fields whose sizes are the Schlafli
   numbers themselves: PSL(2,3) = A_4 (tetrahedron {3,3}), PGL(2,3) = S_4 (cube/octahedron {4,3},{3,4}),
   PSL(2,5) = A_5 = PSL(2,4) = SL(2,4) (icosahedron/dodecahedron {3,5},{5,3}); also PSL(2,2) = S_3 (the triangle).
   Verified here by constructing each group as a permutation group on the projective line P^1(F_q) and identifying
   it: order and parity on 4 points (A_4, S_4), order and parity on 5 points (SL(2,4) = A_5), and for PSL(2,5) the
   conjugation action on its five Klein four-subgroups (faithful, image of order 60 inside A_5).
 * The Fermat primes: 3 = F_0 and 5 = F_1 are the only Fermat primes among Schlafli numbers; PSL(2,17) has order
   2448 and is not a finite rotation group of 3-space (those have orders n, 2n, 12, 24, 60).
 * Constructibility: the icosahedron (0, +-1, +-phi) and dodecahedron (+-1,+-1,+-1), (0, +-1/phi, +-phi) checked
   as regular (equal edges, 5 resp. 3 nearest neighbours); their coordinates lie in Q(sqrt 5), so Euclid XIII's
   constructions use exactly the constructibility of the 3-, 4- and 5-gon.
Usage: python3 gilbreath_fermat_platonic_20260926_groups.py
"""
import itertools, math
import numpy as np


def schlafli():
    sph = [(p, q) for p in range(3, 30) for q in range(3, 30) if (p - 2) * (q - 2) < 4]
    euc = [(p, q) for p in range(3, 30) for q in range(3, 30) if (p - 2) * (q - 2) == 4]
    return sph, euc


class Field:
    def __init__(self, q):
        self.q = q
        if q == 4:
            self.add = lambda x, y: x ^ y
            self.neg = lambda x: x
            self.mul = self._mul4
            self.inv = {1: 1, 2: 3, 3: 2}
        else:
            self.add = lambda x, y: (x + y) % q
            self.neg = lambda x: (-x) % q
            self.mul = lambda x, y: (x * y) % q
            self.inv = {x: pow(x, q - 2, q) for x in range(1, q)}

    @staticmethod
    def _mul4(x, y):
        r = 0
        for i in range(2):
            if (y >> i) & 1:
                r ^= x << i
        if r & 4:
            r ^= 7  # a^2 = a + 1
        return r


def mobius_perm(M, K):
    """permutation of P^1(F_q): points 0..q-1 and q = infinity"""
    q = K.q; a, b, c, d = M
    perm = []
    for x in range(q + 1):
        if x == q:  # infinity
            perm.append(q if c == 0 else K.mul(a, K.inv[c]))
        else:
            den = K.add(K.mul(c, x), d)
            num = K.add(K.mul(a, x), b)
            perm.append(q if den == 0 else K.mul(num, K.inv[den]))
    return tuple(perm)


def closure(gens):
    n = len(gens[0])
    ident = tuple(range(n))
    seen = {ident}; frontier = [ident]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                gh = tuple(h[g[i]] for i in range(n))
                if gh not in seen:
                    seen.add(gh); new.append(gh)
        frontier = new
    return seen


def parity(p):
    n = len(p); seen = [False] * n; sgn = 0
    for i in range(n):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]:
                seen[j] = True; j = p[j]; L += 1
            sgn += L - 1
    return sgn % 2


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(q)))


def inverse(p):
    inv = [0] * len(p)
    for i, x in enumerate(p):
        inv[x] = i
    return tuple(inv)


def order(p):
    ident = tuple(range(len(p))); k = 1; x = p
    while x != ident:
        x = compose(p, x); k += 1
    return k


def psl(q, pgl=False):
    K = Field(q)
    if q == 4:
        gens = [mobius_perm((1, 1, 0, 1), K), mobius_perm((1, 2, 0, 1), K), mobius_perm((0, 1, 1, 0), K)]
    else:
        gens = [mobius_perm((1, 1, 0, 1), K), mobius_perm((0, K.neg(1), 1, 0), K)]
        if pgl:
            ns = next(x for x in range(2, q) if all((y * y) % q != x for y in range(q)))
            gens.append(mobius_perm((1, 0, 0, ns), K))
    return closure(gens)


def klein_subgroups(G):
    ident = tuple(range(len(next(iter(G)))))
    invs = [g for g in G if g != ident and compose(g, g) == ident]
    subs = set()
    for a, b in itertools.combinations(invs, 2):
        if compose(a, b) == compose(b, a):
            subs.add(frozenset([ident, a, b, compose(a, b)]))
    return sorted(subs, key=lambda s: sorted(s)), invs


def main():
    sph, euc = schlafli()
    print("== Schlafli ==")
    print(" (p-2)(q-2) < 4:", sph, "-> %d Platonic solids" % len(sph))
    print(" (p-2)(q-2) = 4:", euc, "-> %d regular plane tilings" % len(euc))
    assert len(sph) == 5 and len(euc) == 3
    print("== projective groups on P^1(F_q) ==")
    for q, pgl, name in ((2, False, 'PSL(2,2)'), (3, False, 'PSL(2,3)'), (3, True, 'PGL(2,3)'), (4, False, 'SL(2,4)'), (5, False, 'PSL(2,5)'), (7, False, 'PSL(2,7)'), (17, False, 'PSL(2,17)')):
        G = psl(q, pgl)
        even = all(parity(g) == 0 for g in G)
        orders = {}
        for g in G:
            o = order(g); orders[o] = orders.get(o, 0) + 1
        print(" %s: order %d on %d points, all even: %s, element orders %s" % (name, len(G), q + 1, even, dict(sorted(orders.items()))))
        if name == 'PSL(2,2)':
            assert len(G) == 6 and not even  # S_3
        if name == 'PSL(2,3)':
            assert len(G) == 12 and even  # A_4: the unique subgroup of order 12 of S_4
        if name == 'PGL(2,3)':
            assert len(G) == 24  # S_4
        if name == 'SL(2,4)':
            assert len(G) == 60 and even  # A_5: the unique subgroup of order 60 of S_5
        if name == 'PSL(2,5)':
            assert len(G) == 60
            subs, invs = klein_subgroups(G)
            print("   PSL(2,5): %d involutions, %d Klein four-subgroups" % (len(invs), len(subs)))
            assert len(subs) == 5
            idx = {s: i for i, s in enumerate(subs)}
            images = set()
            kernel = 0
            for g in G:
                gi = inverse(g)
                img = tuple(idx[frozenset(compose(compose(g, h), gi) for h in s)] for s in subs)
                images.add(img)
                if img == (0, 1, 2, 3, 4):
                    kernel += 1
            assert kernel == 1 and len(images) == 60 and all(parity(p) == 0 for p in images)
            print("   conjugation action on the five Klein four-subgroups: faithful (kernel trivial), image of order 60 of even permutations of 5 objects = A_5")
        if name == 'PSL(2,17)':
            assert len(G) == 2448 and 2448 not in (12, 24, 60)
    print("== constructibility of the solids that use 5 = F_1 ==")
    phi = (1 + 5 ** 0.5) / 2
    ico = [np.array(v) for s1 in (1, -1) for s2 in (1, -1) for v in ((0, s1, s2 * phi), (s1, s2 * phi, 0), (s2 * phi, 0, s1))]
    dod = [np.array(v) for v in itertools.product((1, -1), repeat=3)] + [np.array(v) for s1 in (1, -1) for s2 in (1, -1) for v in ((0, s1 / phi, s2 * phi), (s1 / phi, s2 * phi, 0), (s2 * phi, 0, s1 / phi))]
    for name, V, deg in (('icosahedron', ico, 5), ('dodecahedron', dod, 3)):
        D = np.array([[np.linalg.norm(a - b) for b in V] for a in V])
        e = min(D[D > 1e-9])
        nn = [(np.abs(D[i] - e) < 1e-9).sum() for i in range(len(V))]
        assert len(V) == (12 if deg == 5 else 20) and all(k == deg for k in nn)
        print(" %s: %d vertices, edge %.6f, every vertex has %d nearest neighbours; coordinates in Q(sqrt 5)" % (name, len(V), e, deg))
    print("== the fives ==")
    print(" Fermat primes known: 3, 5, 17, 257, 65537 (F_0..F_4); Schlafli numbers {3,4,5} = {F_0, 2^2, F_1}; PSL fields 3, 4, 5.")
    print(" The only Fermat primes that are Schlafli numbers: 3 and 5. First Fermat prime without a solid: 17 ((17-2)(q-2) >= 15 > 4).")


if __name__ == '__main__':
    main()
