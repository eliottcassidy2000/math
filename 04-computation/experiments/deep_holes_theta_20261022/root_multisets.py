"""(session opus-2026-10-22-S53; THM-4760 (4')) Galois invariance of the root eigenvalue multiset of the N-side shell-1 element
x = (rho_j + theta_j)/(h+2) on the carrier component j and rho_i/(h+2) elsewhere (a necessary condition for (0)).
Roots of each simply-laced component are generated from its Cartan matrix; a root alpha = sum c_i alpha_i has
height sum c_i and theta.alpha = sum c_i (theta.alpha_i)."""
import math
from collections import Counter


def cartan(t, n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    if t == 'A':
        edges = [(i, i + 1) for i in range(n - 1)]
    elif t == 'D':
        edges = [(i, i + 1) for i in range(n - 2)] + [(n - 3, n - 1)]
    elif t == 'E':
        # Bourbaki: 1-3-4-5-...-n chain, 2 attached to 4 (0-based: 0-2-3-4-..., 1 attached to 3)
        edges = [(0, 2), (2, 3), (1, 3)] + [(i, i + 1) for i in range(3, n - 1)]
    for a, b in edges:
        A[a][b] = A[b][a] = -1
    return A


def positive_roots(A):
    n = len(A)
    simple = [tuple(1 if j == i else 0 for j in range(n)) for i in range(n)]
    roots = set(simple)
    frontier = list(simple)
    while frontier:
        new = []
        for r in frontier:
            for i in range(n):
                # alpha_i-string: r + alpha_i is a root iff <r, alpha_i> < 0 or (string) ; use reflection closure
                ip = sum(r[j] * A[j][i] for j in range(n))
                s = list(r); s[i] -= ip
                s = tuple(s)
                if all(c >= 0 for c in s) and any(s) and s not in roots:
                    roots.add(s); new.append(s)
        frontier = new
    return sorted(roots, key=sum)


def comp_data(t, n):
    A = cartan(t, n)
    pos = positive_roots(A)
    theta = pos[-1]
    th_ai = [sum(theta[j] * A[j][i] for j in range(n)) for i in range(n)]
    out = []
    for r in pos:
        ht = sum(r)
        ta = sum(r[i] * th_ai[i] for i in range(n))
        out.append((ht, ta))
    h = sum(theta) + 1
    return out, h


def multiset(components, carrier, m, d):
    c = Counter()
    for k, (t, n) in enumerate(components):
        data, h = comp_data(t, n)
        for ht, ta in data:
            val = ht + (ta if k == carrier else 0)
            c[(d * val) % m] += 1
            c[(-d * val) % m] += 1
    return c


TYPES = {
    'a1^24': [('A', 1)] * 24, 'd4^6': [('D', 4)] * 6, 'd6^4': [('D', 6)] * 4, 'd8^3': [('D', 8)] * 3,
    'd10e7^2': [('D', 10), ('E', 7), ('E', 7)], 'd12^2': [('D', 12)] * 2, 'd16e8': [('D', 16), ('E', 8)],
    'e8^3': [('E', 8)] * 3, 'd24': [('D', 24)],
    'a5^4d4': [('A', 5)] * 4 + [('D', 4)], 'e6^4': [('E', 6)] * 4, 'a15d9': [('A', 15), ('D', 9)],
}
for name, comps in TYPES.items():
    h = comp_data(*comps[0])[1]
    m = h + 2
    units = [d for d in range(1, m) if math.gcd(d, m) == 1]
    seen = set()
    for j, cj in enumerate(comps):
        if cj in seen:
            continue
        seen.add(cj)
        base = multiset(comps, j, m, 1)
        good = [d for d in units if multiset(comps, j, m, d) == base]
        print(f'{name:9s} carrier {cj[0]}{cj[1]:<2d} (m = {m}): root multiset invariant under {len(good)}/{len(units)} units'
              + ('' if len(good) == len(units) else f'  {good}'))


# Level-2 reading: x = (rho + lambda)/(h + 2) with lambda a level-2 weight; carrier = adjoint (theta), others = 2 Lambda_0.
# Galois action on the level-2 primaries of E8 (k + h = 32) and E7 (k + h = 20), compared through root spectra.
def spectrum(t, n, lam, m, d):
    A = cartan(t, n)
    c = Counter()
    for r in positive_roots(A):
        val = sum(r) + sum(lam.get(i, 0) * r[i] for i in range(n))
        c[(d * val) % m] += 1
        c[(-d * val) % m] += 1
    return c


print()
for (t, n, m) in [('E', 8, 32), ('E', 7, 20)]:
    A = cartan(t, n)
    th = positive_roots(A)[-1]
    prim = {'2L0': {}}
    for i in range(n):
        if th[i] == 2:
            prim[f'L{i + 1}'] = {i: 1}
        if th[i] == 1:
            prim[f'L0+L{i + 1}'] = {i: 1}
            for j in range(i, n):
                if th[j] == 1:
                    prim[f'L{i + 1}+L{j + 1}'] = {i: 2} if i == j else {i: 1, j: 1}
    adj = [i for i in range(n) if sum(th[j] * A[j][i] for j in range(n)) != 0][0]
    specs = {k: spectrum(t, n, v, m, 1) for k, v in prim.items()}
    print(f'{t}{n} level 2 (k + h = {m}): primaries {list(prim)}, adjoint L{adj + 1}')
    for d in [d for d in range(2, m // 2) if math.gcd(d, m) == 1]:
        row = []
        for k, v in prim.items():
            img = spectrum(t, n, v, m, d)
            hit = [k2 for k2, s2 in specs.items() if s2 == img]
            row.append(f'{k}->{"/".join(hit) if hit else "?"}')
        print(f'   sigma_{d}: ' + ', '.join(row))
