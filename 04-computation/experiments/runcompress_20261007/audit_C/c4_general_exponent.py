#!/usr/bin/env python3
"""audit_C item 4: decay exponent of the orphan fraction for general residual sources n = 2^K t - 1, as a function of
log2 n ~ K + B, for K = 2 (one child), 3, 5, 9, 17, 33, with a larger sample and a longer range of B (up to 3200 bits).
For K = 2 the orphan event is exactly the failure of the single D = 1 chain, whose exponent is 1/2 by THM-4581 (4)
(sketch): this calibrates the method. Usage: python3 c4_general_exponent.py N seed K1,K2,...
"""
import random, sys, math
exec(open('c4_general_sources.py').read().split('N = int(sys.argv[1])')[0])   # reuse ot(), v2()
N = int(sys.argv[1]); seed = int(sys.argv[2]); Ks = [int(x) for x in sys.argv[3].split(',')]
rnd = random.Random(seed)
Bs = (25, 50, 100, 200, 400, 800, 1600, 3200)
for K in Ks:
    pts = []
    for B in Bs:
        orph = 0; tot = 0
        while tot < N:
            t = rnd.getrandbits(B) | (1 << (B - 1)) | 1
            if v2(3 ** K * t - 1) != 1:
                continue
            n = (t << K) - 1
            on, tn = ot(n)
            found = False
            for D in range(1, K):
                oh, th = ot(((n + 1) >> D) - 1)
                if oh == on and tn - th == D:
                    found = True; break
            tot += 1; orph += not found
        pts.append((K + B, orph, N))
        print(f"K={K:2d} B={B:5d}: orphans {orph}/{N} = {orph/N:.4f}, x sqrt(K+B) = {orph/N*math.sqrt(K+B):.3f}", flush=True)
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
    for lab, sub in (("all B", pts), ("B >= 200", [p for p in pts if p[0] - K >= 200])):
        prof = [(0.2 + 0.005 * i, ll(0.2 + 0.005 * i, sub)) for i in range(201)]
        best = max(prof, key=lambda z: z[1]); ci = [a for a, l in prof if l >= best[1] - 1.92]
        print(f"K={K} ({lab}): alpha_hat = {best[0]:.3f}, 95% profile interval [{min(ci):.3f}, {max(ci):.3f}], "
              f"log-lik drop at 1/2 = {best[1] - ll(0.5, sub):.2f}", flush=True)
