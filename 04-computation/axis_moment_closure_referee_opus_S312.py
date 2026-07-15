# opus-2026-07-15-S312 -- F6 referee: the axis walk second-moment closure.
# Claims (N = C(n,2), x = sum_v d_v^2, d = 2s - (n-1)):
#   F5:  E[Delta x | T]   = 8 - 4x/N                       (re-verified)
#   F6a: E[Delta^2 | T]   = 16(n-4)x/N + 64                (NEW; exact, per T)
#   F6b: sum_pairs (d_u - d_v)^2 = n*x                     (the quadratic identity)
#   F6c: stationary Var(x) = 2n(n-1)(n-2)                  (uniform tournaments)
#   F6d: NO third-moment closure: E[Delta^3 | T] is NOT a function of x alone
# Exhaustive n = 4, 5, 6.
from collections import defaultdict
from fractions import Fraction

for n in (4, 5, 6):
    N = n*(n-1)//2
    pairs = [(u, v) for u in range(n) for v in range(u+1, n)]
    ok5 = ok6a = ok6b = True
    xs = []
    third_by_x = defaultdict(set)
    for X in range(1 << N):
        s = [0]*n
        for i, (u, v) in enumerate(pairs):
            if (X >> i) & 1: s[u] += 1
            else: s[v] += 1
        d = [2*sv - (n-1) for sv in s]
        x = sum(dv*dv for dv in d)
        xs.append(x)
        # per-arc flip deltas
        deltas = []
        for i, (u, v) in enumerate(pairs):
            if (X >> i) & 1: w, l = u, v
            else: w, l = v, u
            deltas.append(4*(d[l] - d[w]) + 8)
        E1 = Fraction(sum(deltas), N)
        E2 = Fraction(sum(dd*dd for dd in deltas), N)
        E3 = Fraction(sum(dd**3 for dd in deltas), N)
        if E1 != 8 - Fraction(4*x, N): ok5 = False
        if E2 != Fraction(16*(n-4)*x, N) + 64: ok6a = False
        q = sum((d[u]-d[v])**2 for (u, v) in pairs)
        if q != n*x: ok6b = False
        third_by_x[x].add(E3)
    mean = Fraction(sum(xs), len(xs))
    var = Fraction(sum(v*v for v in xs), len(xs)) - mean*mean
    split_x = [x for x, s3 in third_by_x.items() if len(s3) > 1]
    print(f"n={n}: F5 drift {ok5} | F6a second-moment {ok6a} | F6b quadratic {ok6b}")
    print(f"   E[x] = {mean} (predict {n*(n-1)}) | Var(x) = {var} "
          f"(predict {2*n*(n-1)*(n-2)}) -> {var == 2*n*(n-1)*(n-2)}")
    print(f"   F6d third-moment split x-values: {len(split_x)} "
          f"(closure fails at {sorted(split_x)[:4]}{'...' if len(split_x)>4 else ''})"
          if split_x else "   F6d: E[Delta^3] constant per x at this n (no witness)")
