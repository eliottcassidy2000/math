# list m with supersingular density >= 0.3 for every family, |m| <= M (vectorised by residues)
import sys, json, numpy as np
sys.path.insert(0, '.')
from sporadic_cm_detector import FAMILIES, primes_between
M = int(sys.argv[1]); ps = primes_between(7, 500)
ms = np.array([m for m in range(-M, M+1) if m != 0], dtype=np.int64)
out = {}
for name, gen in FAMILIES.items():
    u = gen(); hits = np.zeros(len(ms)); tot = np.zeros(len(ms))
    for p in ps:
        # zero set of sum_k u_k r^{-k} over residues r: evaluate polynomial sum_k u_k y^k at y = r^{-1}
        coeffs = np.array([x % p for x in u[:p]], dtype=np.int64)
        y = np.arange(p, dtype=np.int64); val = np.zeros(p, dtype=np.int64)
        for c in coeffs[::-1]:
            val = (val * y + c) % p
        inv = np.zeros(p, dtype=np.int64)
        for r in range(1, p): inv[r] = pow(r, -1, p)
        zero_r = (val[inv] == 0); zero_r[0] = False
        res = np.mod(ms, p); ok = res != 0
        tot += ok; hits += ok & zero_r[res]
    d = hits / tot
    sel = ms[d >= 0.3]
    out[name] = {int(m): round(float(x), 3) for m, x in zip(ms[d >= 0.3], d[d >= 0.3])}
    print(name, out[name], flush=True)
json.dump(out, open(f'highdensity_M{M}.json', 'w'), indent=1)
