#!/usr/bin/env python3
"""audit_C item 4: general residual sources n = 2^K t - 1 (t odd, B bits, first reset letter 2: v2(3^K t - 1) = 1).
Orphan: no child h_D = (n+1)/2^D - 1, 1 <= D <= K-1, with equal odd count and Terras difference D.
Independent rerun (own orbit code with 16-step residue tables; different seeds) with N sources per row, plus a binomial MLE of
the decay exponent in log2 n ~ K + B, with profile-likelihood intervals.
Usage: python3 c4_general_sources.py N seed
"""
import random, sys, math

BB = 16
MASK = (1 << BB) - 1
ODD = [0] * (1 << BB)
VAL = [0] * (1 << BB)
for r in range(1 << BB):
    x, o = r, 0
    for _ in range(BB):
        if x & 1:
            x = (3 * x + 1) >> 1
            o += 1
        else:
            x >>= 1
    ODD[r], VAL[r] = o, x
P3 = [3 ** i for i in range(BB + 1)]
LIM = 1 << (BB + 1)


def ot(x):
    odd = ter = 0
    while x >= LIM:
        r = x & MASK
        o = ODD[r]
        x = P3[o] * (x >> BB) + VAL[r]
        odd += o
        ter += BB
    while x != 1:
        if x & 1:
            x = (3 * x + 1) >> 1
            odd += 1
        else:
            x >>= 1
        ter += 1
    return odd, ter


def v2(x):
    return (x & -x).bit_length() - 1


N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 11
rnd = random.Random(seed)
Bs = (25, 50, 100, 200, 400, 800)
results = {}
for K in (9, 17):
    for B in Bs:
        orph = d12 = 0
        tot = 0
        while tot < N:
            t = rnd.getrandbits(B) | (1 << (B - 1)) | 1
            if v2(3 ** K * t - 1) != 1:
                continue
            n = (t << K) - 1
            on, tn = ot(n)
            found = False
            f12 = False
            for D in range(1, K):
                h = ((n + 1) >> D) - 1
                oh, th = ot(h)
                if oh == on and tn - th == D:
                    found = True
                    if D <= 2:
                        f12 = True
                    break
            tot += 1
            orph += not found
            d12 += f12
        results[(K, B)] = (orph, N)
        print(f"K={K:2d} B={B:4d}: orphan fraction {orph/N:.3f} ({orph}/{N}), x sqrt(K+B) = {orph/N*math.sqrt(K+B):.2f}; "
              f"reset pair D<=2 merges {d12/N:.3f}", flush=True)

# binomial MLE of p = C (K+B)^-alpha per K
for K in (9, 17):
    pts = [(K + B, results[(K, B)][0], results[(K, B)][1]) for B in Bs]

    def ll(alpha):
        # profile over C (closed form not available because of the clip; golden search on log C)
        def f(lc):
            s = 0.0
            for x, k, n in pts:
                p = min(max(math.exp(lc) * x ** (-alpha), 1e-12), 1 - 1e-12)
                s += k * math.log(p) + (n - k) * math.log(1 - p)
            return s
        a, b = -10.0, 10.0
        g = (math.sqrt(5) - 1) / 2
        c, d = b - g * (b - a), a + g * (b - a)
        for _ in range(100):
            if f(c) > f(d):
                b = d
            else:
                a = c
            c, d = b - g * (b - a), a + g * (b - a)
        return f((a + b) / 2)
    prof = [(0.2 + 0.005 * i, ll(0.2 + 0.005 * i)) for i in range(161)]
    best = max(prof, key=lambda z: z[1])
    ci = [a for a, l in prof if l >= best[1] - 1.92]
    l_half = ll(0.5)
    print(f"K={K}: alpha_hat = {best[0]:.3f}, 95% profile interval [{min(ci):.3f}, {max(ci):.3f}]; "
          f"log-lik drop at alpha = 1/2: {best[1] - l_half:.2f}")
