# Identify the arithmetic term in TT5(P_p): least-squares polynomial main term (degree 5 in p), then the residual's
# size scaling and its correlation with CM data for primes p = 3 (mod 4).
import re, math
import numpy as np
rows = [(int(m.group(1)), int(m.group(2))) for m in (re.match(r'p=(\d+) TT5=(\d+)', l) for l in open('paley_tt5_data.out')) if m]
P = np.array([r[0] for r in rows], dtype=float); T = np.array([r[1] for r in rows], dtype=object)
# exact rational fit would be ideal; use float on TT5 / p^5 scale
y = np.array([float(t) for t in T])
V = np.vander(P, 6)
coef, *_ = np.linalg.lstsq(V, y, rcond=None)
res = y - V @ coef
print("leading coefficients (p^5..p^0):", [f"{c:.6g}" for c in coef])
for lo, hi in ((19, 300), (300, 800), (800, 1500)):
    sel = (P >= lo) & (P < hi)
    print(f"  p in [{lo},{hi}): rms residual {np.sqrt(np.mean(res[sel] ** 2)):.4g}; rms/p^2 {np.sqrt(np.mean((res[sel] / P[sel] ** 2) ** 2)):.4g}; rms/p^2.5 {np.sqrt(np.mean((res[sel] / P[sel] ** 2.5) ** 2)):.4g}")
def rep(p, a, b, c):
    """some x with p = a x^2 + b x y + c y^2 (returns |x| list over representations) or None"""
    out = []
    lim = int(math.isqrt(4 * c * p // max(1, 4 * a * c - b * b))) + 2
    for x in range(-lim, lim + 1):
        # solve c y^2 + b x y + a x^2 - p = 0 for integer y
        A_, B_, C_ = c, b * x, a * x * x - p
        disc = B_ * B_ - 4 * A_ * C_
        if disc < 0: continue
        s = math.isqrt(disc)
        if s * s != disc: continue
        for yv in ((-B_ + s), (-B_ - s)):
            if yv % (2 * A_) == 0: out.append((x, yv // (2 * A_)))
    return out
cands = {}
for name, (a, b, c) in {'x^2+2y^2': (1, 0, 2), 'x^2+7y^2': (1, 0, 7), 'x^2+xy+2y^2': (1, 1, 2), 'x^2+3y^2': (1, 0, 3),
                        'x^2+xy+y^2': (1, 1, 1), 'x^2+11y^2': (1, 0, 11), 'x^2+xy+3y^2': (1, 1, 3)}.items():
    f1 = []; f2 = []
    for p in P.astype(int):
        r = rep(p, a, b, c)
        if not r: f1.append(0.0); f2.append(0.0); continue
        x, yy = r[0]
        f1.append(float(x * x)); f2.append(float(x * x - c * yy * yy * (1 if b == 0 else 1)))
    cands[name + ' (x^2)'] = np.array(f1); cands[name + ' (x^2 - c y^2)'] = np.array(f2)
z = res / P ** 2           # normalise by p^2 (if the arithmetic term is p^2 * (weight-2-like number of size sqrt p) ...)
for name, f in cands.items():
    if np.std(f) == 0: continue
    for scale_name, g in (('res', res), ('res/p', res / P), ('res/p^2', res / P ** 2)):
        cc = np.corrcoef(g, f)[0, 1]
        if abs(cc) > 0.5: print(f"  corr({scale_name}, {name}) = {cc:+.3f}")
print("residual / p^2 at the first 12 primes:", [f"{v:+.2f}" for v in (res / P ** 2)[:12]])
print("p mod 8 / mod 3 / mod 7 of those:", [(int(p) % 8, int(p) % 3, int(p) % 7) for p in P[:12]])
