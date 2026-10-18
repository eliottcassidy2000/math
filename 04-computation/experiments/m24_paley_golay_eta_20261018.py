#!/usr/bin/env python3
"""From the Paley tournament P_23 to the Golay code, M24, and Mason's eta-products (S43).

1. The quadratic residues Q mod 23 (the out-neighbourhood of 0 in the Paley tournament P_23) generate the binary
   extended quadratic-residue code of length 24 on P^1(F_23) = {0..22, oo}; we check that it is the Golay code
   (dimension 12, weight enumerator 1 + 759 y^8 + 2576 y^12 + 759 y^16 + y^24).
2. M24 = <alpha, beta, gamma, delta> on P^1(F_23) (Conway's generators: x+1, 2x, -1/x, and delta: x -> x^3/9 on
   residues, 9 x^3 on non-residues, fixing 0 and oo); we check that every generator preserves the code, and sample
   (approximately uniform) group elements by long random words to census the cycle shapes against the 26 classes,
   comparing frequencies with 1/|C(g)|.
3. For each of the 26 cycle shapes prod a^(n_a) (sum a n_a = 24) the eta-product eta_g = prod eta(a tau)^(n_a)
   (weight = #cycles/2) is expanded; multiplicativity a(mn) = a(m)a(n) (coprime m, n) and the Hecke relation at primes
   a(p^2) = a(p)^2 - chi(p) p^(k-1) are tested, the latter with the nebentypus fitted as +-1.
Usage: python m24_paley_golay_eta_20261018.py
Session opus-2026-10-17-S43 (date 2026-10-18).
"""
import sys, os, json, random, math
from fractions import Fraction

P = 23
INF = 23
QR = sorted({(x * x) % P for x in range(1, P)})
NQR = sorted(set(range(1, P)) - set(QR))

SHAPES = {  # the 26 conjugacy classes of M24 and their cycle shapes; centralizer orders from the ATLAS
    '1A': ({1: 24}, 244823040), '2A': ({1: 8, 2: 8}, 21504), '2B': ({2: 12}, 7680), '3A': ({1: 6, 3: 6}, 1080),
    '3B': ({3: 8}, 504), '4A': ({2: 4, 4: 4}, 384), '4B': ({1: 4, 2: 2, 4: 4}, 128), '4C': ({4: 6}, 96),
    '5A': ({1: 4, 5: 4}, 60), '6A': ({1: 2, 2: 2, 3: 2, 6: 2}, 24), '6B': ({6: 4}, 24), '7AB': ({1: 3, 7: 3}, 42),
    '8A': ({1: 2, 2: 1, 4: 1, 8: 2}, 16), '10A': ({2: 2, 10: 2}, 20), '11A': ({1: 2, 11: 2}, 11),
    '12A': ({2: 1, 4: 1, 6: 1, 12: 1}, 12), '12B': ({12: 2}, 12), '14AB': ({1: 1, 2: 1, 7: 1, 14: 1}, 14),
    '15AB': ({1: 1, 3: 1, 5: 1, 15: 1}, 15), '21AB': ({3: 1, 21: 1}, 21), '23AB': ({1: 1, 23: 1}, 23),
}
PAIRED = {'7AB', '14AB', '15AB', '21AB', '23AB'}   # two classes (A, B) with the same cycle shape


def inv(x):
    return pow(x, P - 2, P)


def perm_from(fun):
    return tuple(fun(x) for x in range(24))


def alpha(x):
    return INF if x == INF else (x + 1) % P


def beta(x):
    return INF if x == INF else (2 * x) % P


def gamma(x):
    if x == INF:
        return 0
    if x == 0:
        return INF
    return (-inv(x)) % P


def delta(x):
    if x == INF or x == 0:
        return x
    if x in QR:
        return (pow(x, 3, P) * inv(9)) % P
    return (9 * pow(x, 3, P)) % P


def compose(p, q):   # p after q
    return tuple(p[q[i]] for i in range(24))


def cycle_shape(p):
    seen = [False] * 24
    sh = {}
    for i in range(24):
        if not seen[i]:
            L = 0
            j = i
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            sh[L] = sh.get(L, 0) + 1
    return sh


def golay():
    gens = []
    for a in range(P):
        v = [0] * 24
        for r in [0] + QR:
            v[(a + r) % P] = 1
        gens.append(v)
    gens.append([1] * 24)
    # row reduce over F2 (bitmasks)
    rows = [sum(b << i for i, b in enumerate(v)) for v in gens]
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r:
            basis.append(r)
    # canonical: full reduction not needed; enumerate span
    code = {0}
    for b in basis:
        code |= {c ^ b for c in code}
    return basis, code


def weight_enumerator(code):
    we = {}
    for c in code:
        w = bin(c).count('1')
        we[w] = we.get(w, 0) + 1
    return dict(sorted(we.items()))


def preserves(perm, code):
    for c in code:
        img = 0
        for i in range(24):
            if (c >> i) & 1:
                img |= 1 << perm[i]
        if img not in code:
            return False
    return True


def shape_key(sh):
    for name, (s, cen) in SHAPES.items():
        if s == sh:
            return name
    return 'OTHER:' + str(sorted(sh.items()))


def eta_product(shape, N):
    """coefficients b[0..N] with eta_g = q^(sum a n_a / 24) * sum b[n] q^n."""
    b = [0] * (N + 1)
    b[0] = 1
    for a, na in shape.items():
        for _ in range(na):
            for n in range(1, N // a + 1):
                step = a * n
                for m in range(N, step - 1, -1):
                    b[m] -= b[m - step]
    return b


def main():
    out = {}
    basis, code = golay()
    we = weight_enumerator(code)
    print('extended QR code from Paley P_23 residues: dimension %d, weight enumerator %s' % (len(basis), we))
    out['golay'] = dict(dimension=len(basis), weight_enumerator=we)
    gens = {'alpha': perm_from(alpha), 'beta': perm_from(beta), 'gamma': perm_from(gamma), 'delta': perm_from(delta)}
    pres = {k: preserves(v, code) for k, v in gens.items()}
    print('generators preserve the code:', pres)
    out['generators_preserve'] = pres
    # Paley automorphisms x -> a x + b (a in QR) as M24 elements
    aut = {}
    for a in QR:
        for bb in range(P):
            pm = perm_from(lambda x, a=a, bb=bb: INF if x == INF else (a * x + bb) % P)
            k = shape_key(cycle_shape(pm))
            aut[k] = aut.get(k, 0) + 1
    print('Aut(P_23) = {x -> a x + b, a in QR} (order 253) inside M24, by class:', aut)
    out['paley_aut_classes'] = aut
    # random sampling
    random.seed(20261018)
    G = list(gens.values())
    cur = tuple(range(24))
    for _ in range(2000):
        cur = compose(random.choice(G), cur)
    counts = {}
    S = 400000
    for t in range(S):
        for _ in range(6):
            cur = compose(random.choice(G), cur)
        k = shape_key(cycle_shape(cur))
        counts[k] = counts.get(k, 0) + 1
    rows = []
    print('cycle-shape census over %d sampled elements (random walk, 6 steps apart):' % S)
    for name, (sh, cen) in SHAPES.items():
        mult = 2 if name in PAIRED else 1
        exp = mult / cen
        obs = counts.get(name, 0) / S
        rows.append(dict(cls=name, shape=sh, centralizer=cen, expected=exp, observed=obs))
        print('   %-5s %-28s expected %.5f observed %.5f' % (name, ' '.join('%d^%d' % (a, n) for a, n in sorted(sh.items())), exp, obs))
    other = {k: v for k, v in counts.items() if k.startswith('OTHER')}
    print('   shapes outside the 26 classes:', other)
    out['census'] = dict(samples=S, rows=rows, other=other)
    # eta products
    N = 600
    etas = {}
    primes = [p for p in range(2, N + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    print('Mason eta-products eta_g (weight = #cycles/2):')
    for name, (sh, cen) in SHAPES.items():
        k2 = sum(sh.values())          # 2 * weight
        level = math.lcm(*sh.keys()) * min(sh.keys())   # heuristic level N_g = (largest * smallest cycle) = ord(g) * h
        b = eta_product(sh, N)
        off = sum(a * n for a, n in sh.items()) // 24   # = 1
        a = [0] * (N + 1)
        for n in range(N + 1 - off):
            a[n + off] = b[n]
        mult_ok = all(a[m * n] == a[m] * a[n] for m in range(2, 40) for n in range(2, 40)
                      if math.gcd(m, n) == 1 and m * n <= N)
        hecke = []
        wt = k2 // 2 if k2 % 2 == 0 else None
        if wt is not None:
            for p in primes:
                if p * p > N:
                    break
                if level % p == 0:
                    continue
                d = a[p] ** 2 - a[p * p]
                hecke.append(Fraction(d, p ** (wt - 1)))
        chis = sorted(set(hecke))
        etas[name] = dict(shape=sh, weight=k2 / 2, level_guess=level, coeffs=a[:61], multiplicative=mult_ok,
                          hecke_chi_values=[str(c) for c in chis])
        print('   %-5s weight %-4s level %-4d multiplicative %-5s Hecke chi(p) values %s; a(1..12) = %s'
              % (name, k2 / 2, level, mult_ok, [str(c) for c in chis], a[1:13]))
    out['eta_products'] = etas
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'm24_paley_golay_eta_20261018.json'), 'w'),
              indent=1, default=str)


if __name__ == '__main__':
    main()
