#!/usr/bin/env python3
"""audit_C item 4: binomial MLE of the decay exponent on the SESSION's own general-source table (orphan_law_general.out,
300 sources per row), p = C (K+B)^-alpha, profile-likelihood 95% interval."""
import math, re
pts = {9: [], 17: []}
for line in open('../orphan_law_general.out'):
    m = re.match(r'K=\s*(\d+) B=\s*(\d+): orphan fraction ([0-9.]+)', line)
    if m:
        K, B, f = int(m.group(1)), int(m.group(2)), float(m.group(3))
        pts[K].append((K + B, round(f * 300), 300))
def ll(alpha, sub):
    def f(lc):
        s = 0.0
        for x, k, n in sub:
            p = min(max(math.exp(lc) * x ** (-alpha), 1e-12), 1 - 1e-12)
            s += k * math.log(p) + (n - k) * math.log(1 - p)
        return s
    a, b = -10.0, 10.0; g = (math.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    for _ in range(100):
        if f(c) > f(d): b = d
        else: a = c
        c, d = b - g * (b - a), a + g * (b - a)
    return f((a + b) / 2)
for K, sub in pts.items():
    prof = [(0.2 + 0.005 * i, ll(0.2 + 0.005 * i, sub)) for i in range(161)]
    best = max(prof, key=lambda z: z[1]); ci = [a for a, l in prof if l >= best[1] - 1.92]
    print(f"session data K={K}: counts {[k for _, k, _ in sub]}/300; alpha_hat = {best[0]:.3f}, 95% interval [{min(ci):.3f}, {max(ci):.3f}], "
          f"log-lik drop at 1/2 = {best[1] - ll(0.5, sub):.2f}")
