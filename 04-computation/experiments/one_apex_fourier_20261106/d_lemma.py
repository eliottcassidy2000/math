# THM-4887 checks (opus-2026-11-06-S83):
# (i) strip lemma with the worst-case shift D: window W = C: Z* <= 2C^2/P + 3C/beta + 5;
#     window W = 2C: Z* <= 4C^2/P + 3C/beta + 7 (beta = reduced Dirichlet denominator of h/h' with parameter B = P/C);
# (ii) the small-denominator count: #{pairs with beta(h,h') = beta} <= P (beta+1)(4C/beta + 1), and
#      sum over h != h' of 1/beta <= P (4C ln B + 10.58 C + B + 1 + ln B);
# (iii) same h, l != l': |g(v) - g(v')| > C;  (iv) E_n(Delta) for Delta in {0, C} against PC + sum Z* and (D1).
import numpy as np, math
from fractions import Fraction
from collections import Counter
def g(v):
    acc = 0
    while v > 0:
        if v & 1: v >>= 1
        else:
            w = v >> 1; acc += w + (w & 1); v = w
    return 1 + acc
def beta_of(h, hp, B):
    for b in range(1, B + 1):
        x = Fraction(b * h, hp)
        a = round(x)
        if abs(x - a) < Fraction(1, B):
            return Fraction(a, b).denominator
    raise RuntimeError
def maxwin(vals, W):
    j = np.searchsorted(vals, vals + 2 * W, side='right')
    return int((j - np.arange(len(vals))).max())
def run(n, m):
    C = 1 << m; k = n - m; P = 1 << (k - 1); B = P // C
    assert P >= 8 * C
    c = np.arange(C, dtype=np.int64)
    viol1 = viol2 = 0; r1 = r2 = 0.0; S1 = S2 = 0; invb = 0.0; nb = Counter()
    for h in range(P, 2 * P):
        hc = h * c
        for hp in range(P, 2 * P):
            if hp == h: continue
            vals = np.sort((hc[:, None] - hp * c[None, :]).ravel())
            z1 = maxwin(vals, C); z2 = maxwin(vals, 2 * C)
            b = beta_of(h, hp, B); nb[b] += 1; invb += 1.0 / b
            b1 = 2 * C * C / P + 3 * C / b + 5; b2 = 4 * C * C / P + 3 * C / b + 7
            viol1 += z1 > b1; viol2 += z2 > b2; r1 = max(r1, z1 / b1); r2 = max(r2, z2 / b2)
            S1 += z1; S2 += z2
    cntviol = sum(1 for b, k_ in nb.items() if k_ > P * (b + 1) * (4 * C / b + 1))
    invb_bound = P * (4 * C * math.log(B) + 10.58 * C + B + 1 + math.log(B))
    # same-h separation and energies
    v = np.arange(1 << (n - 1), 1 << n, dtype=np.int64)
    G = np.array([g(int(x)) for x in v], dtype=np.int64)
    sep = min(int(np.diff(np.sort(G[i * C:(i + 1) * C])).min()) for i in range(P))
    x = np.sort(G)
    E0 = int((np.unique(x, return_counts=True)[1].astype(np.int64) ** 2).sum())
    EC = int((np.searchsorted(x, x + C, side='right') - np.searchsorted(x, x - C, side='left')).sum())
    lB = math.log(B)
    D1_0 = P * C + 8 * P * P + P * C * C * (33.74 + 12 * lB) + 3 * P * C * (1 + lB)
    D1_C = P * C + 10 * P * P + P * C * C * (35.74 + 12 * lB) + 3 * P * C * (1 + lB)
    print(f"n={n} m={m} P={P} C={C} B={B}: strip lemma W=C viol {viol1} (max ratio {r1:.3f}), W=2C viol {viol2} "
          f"(max ratio {r2:.3f}); beta-count violations {cntviol}; sum 1/beta = {invb:.1f} <= {invb_bound:.1f}; "
          f"min same-h gap {sep} > C = {C}; E_n(0) = {E0} <= PC + sum Z*(W=C) = {P*C + S1} <= (D1) {D1_0:.0f}; "
          f"E_n(C) = {EC} <= PC + sum Z*(W=2C) = {P*C + S2} <= (D1) {D1_C:.0f}", flush=True)
for (n, m) in [(10, 2), (11, 2), (12, 3), (13, 3)]:
    run(n, m)
