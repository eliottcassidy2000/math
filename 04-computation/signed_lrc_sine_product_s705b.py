"""SIGNED LRC homometry -- the SINE-PRODUCT closed form + Hamming decomposition.
monad-explorer-2026-06-06-S705.

CLAIM (to verify):  for AP_n, C=2n-1, cut P (=runners colored +), N=complement,
   |f_hat_eps(t)|^2 = |F(t)|^2 - 4*sigma_P(t)*sigma_N(t),
   sigma_X(t) = sum_{i in X} sin(2*pi*t*i/C),  F(t)=sum_{i=1}^{n-1} omega^{t i}.
=>  cuts P,P' COLLIDE  <=>  sigma_P(t)*sigma_N(t) == sigma_P'(t)*sigma_N'(t)  for all t=1..n-1.

Also: decompose deficiency by HAMMING DISTANCE of colliding sign-vectors (how many runners
differ), to separate single-flip (THM-413 order-3) from genuine multi-flip homometries.
"""
from itertools import product, combinations
import cmath, math, sys
from collections import Counter


def folded_spectrum(V, eps, C):
    out = []
    n1 = len(V)
    for i in range(n1):
        for j in range(i + 1, n1):
            f = (eps[i] * V[i] - eps[j] * V[j]) % C
            out.append(min(f, C - f))
    return tuple(sorted(out))


def dft_power_vec(eps, V, C):
    pts = [(eps[i] * V[i]) % C for i in range(len(V))]
    return tuple(round((abs(sum(cmath.exp(2j*math.pi*t*p/C) for p in pts)))**2, 6) for t in range(C))


def sine_product_vec(eps, V, C):
    """[sigma_P(t)*sigma_N(t) for t in 1..(C-1)//2], the cut-dependent part."""
    P = [V[i] for i in range(len(V)) if eps[i] == 1]
    N = [V[i] for i in range(len(V)) if eps[i] == -1]
    out = []
    for t in range(1, (C - 1)//2 + 1):
        sp = sum(math.sin(2*math.pi*t*i/C) for i in P)
        sn = sum(math.sin(2*math.pi*t*i/C) for i in N)
        out.append(round(sp*sn, 6))
    return tuple(out)


def Fpower(V, C):
    return [round((abs(sum(cmath.exp(2j*math.pi*t*i/C) for i in V)))**2, 6) for t in range(C)]


print("=" * 92)
print("PART A: verify  |f_hat|^2 = |F|^2 - 4 sigma_P sigma_N  (exhaustive over cuts)")
print("=" * 92)
for n in range(4, 12):
    V = list(range(1, n)); C = 2*n-1
    Fp = Fpower(V, C)
    ok = True
    for bits in product([0,1], repeat=n-2):
        eps = [1] + [1 if b else -1 for b in bits]
        P = [V[i] for i in range(len(V)) if eps[i]==1]
        N = [V[i] for i in range(len(V)) if eps[i]==-1]
        lhs = dft_power_vec(eps, V, C)
        for t in range(C):
            sp = sum(math.sin(2*math.pi*t*i/C) for i in P)
            sn = sum(math.sin(2*math.pi*t*i/C) for i in N)
            rhs = Fp[t] - 4*sp*sn
            if abs(lhs[t]-rhs) > 1e-6:
                ok = False; break
        if not ok: break
    print(f"  n={n:>2} C={C:>3}: formula holds for all {2**(n-2)} cuts: {ok}", flush=True)

print("\n" + "=" * 92)
print("PART B: collision <=> equal sine-product vector; deficiency by Hamming distance")
print("=" * 92)
print(f"{'n':>3} {'C':>4} {'C-fact':>9} {'defic':>6} {'sine=fold?':>10}   Hamming-dist histogram of colliding pairs")
def factor(m):
    f=[]; d=2
    while d*d<=m:
        while m%d==0: f.append(d); m//=d
        d+=1
    if m>1: f.append(m)
    return f
for n in range(5, 21):
    V = list(range(1, n)); C = 2*n-1
    fold_map = {}
    sine_map = {}
    allcuts = []
    for bits in product([0,1], repeat=n-2):
        eps = tuple([1] + [1 if b else -1 for b in bits])
        allcuts.append(eps)
        fold_map.setdefault(folded_spectrum(V, list(eps), C), []).append(eps)
        sine_map.setdefault(sine_product_vec(list(eps), V, C), []).append(eps)
    # do fold and sine give same partition?
    def part(m): return frozenset(frozenset(g) for g in m.values())
    same = (part(fold_map) == part(sine_map)) if n <= 13 else "(skip)"
    defic = len(allcuts) - len(fold_map)
    # Hamming histogram: within each colliding group, all pairwise hamming dists
    ham = Counter()
    for g in fold_map.values():
        if len(g) > 1:
            for a, b in combinations(g, 2):
                d = sum(1 for i in range(len(V)) if a[i] != b[i])
                d = min(d, len(V)-d)  # global-swap canonical
                ham[d] += 1
    fac = factor(C)
    facstr = "x".join(map(str,fac)) if len(fac)>1 else f"{C}P"
    print(f"{n:>3} {C:>4} {facstr:>9} {defic:>6} {str(same):>10}   {dict(sorted(ham.items()))}", flush=True)
