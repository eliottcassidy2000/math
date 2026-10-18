import sympy as sp
from b2_degree192 import place_values, leg
cnt = {'D': 0, 'R': 0, 'DR': 0}; rows = 0; bp = {}
for p in sp.primerange(29, 3000):
    if p % 15 != 1:
        continue
    js = sp.sqrt_mod(-11 % p, p, all_roots=True)
    ws = [x for x in range(p) if (x ** 3 + 2 * x * x - 4 * x - 4) % p == 0]
    zs = [x for x in range(2, p) if pow(x, 15, p) == 1 and pow(x, 3, p) != 1 and pow(x, 5, p) != 1]
    if not (js and ws and zs):
        continue
    for z0 in zs:
        for j0 in js:
            for w0 in ws:
                N, D, R = place_values(p, z0, j0, w0)
                rows += 1
                for name, val in (('D', D), ('R', R), ('DR', D * R)):
                    if leg(val, p) == -1:
                        cnt[name] += 1
                        bp.setdefault(name, set()).add(p)
print('rows', rows, 'non-residue counts', cnt)
print('primes with a D non-residue:', sorted(bp['D']))
