#!/usr/bin/env python3
"""collatz_procgen_20260922_endgame_diophantine.py -- Diophantine inputs of the Q2 endgame and the renewal
structure of the canonical escape Psi.

  K  LTE in exact 3-adic form: v_3(2^n w - 1) = [n = n0 mod 2] (1 + v_3((n-n0)/2 - alpha_w)) with the 3-adic
     discrete logarithm alpha_w = -log(w')/log(4), w' = 2^n0 w (checked); depth statistics of alpha_w.
  L  Yu (Forum Math. 19 (2007), Thm quoted as A.2 in arXiv:2107.00971): explicit bound for v_3(2^n w - p).
  M  double exhaustion (Senge-Straus 1973 / Stewart 1980): solutions of 2^(a+1) u - 1 = 3^k w' and
     2^a u - 1 = 3^k w' with small u, w' (a landing 2^a u that is again an exhausting hostile point).
  N  the digit-exhausted family m = (1 + 3^k w)/2 (w small): Psi-orbits far beyond 10^18 (exact big integers).
  O  renewal structure of Psi under Haar: drift per digit, dim of the Psi-exceptional set (pressure),
     Lundberg exponent (random-model size of the worst excursion), exact counts of Psi-alive classes.
"""
import sys, math, random
from functools import lru_cache

LN2, LN3 = math.log(2.0), math.log(3.0)
def K0(s): return (3 ** (s + 1)).bit_length() - 1
def Kstar(s): return 0 if s == 0 else (2 if s == 1 else K0(s))
def v3int(n):
    if n == 0: return 10 ** 9
    v = 0
    while n % 3 == 0: n //= 3; v += 1
    return v

def dlog4(wp, N):
    """alpha with 4^alpha * wp = 1 mod 3^(N+1), alpha mod 3^N  (wp = 1 mod 3)."""
    mod = 3 ** (N + 1); M = 0; u = wp % mod
    for t in range(N):
        g = pow(4, 3 ** t, mod)
        for d in range(3):
            if (u * pow(g, d, mod) - 1) % 3 ** (t + 2) == 0:
                M += d * 3 ** t; u = u * pow(g, d, mod) % mod; break
        else: raise AssertionError("dlog")
    return M

def section_K():
    print("=" * 100)
    print("K. LTE in 3-adic form and the digits of alpha_w = -log(w')/log 4")
    N = 40
    nchk = 0; rows = []
    for w in range(1, 400):
        if w % 3 == 0: continue
        n0 = 0 if w % 3 == 1 else 1
        wp = (2 ** n0) * w
        al = dlog4(wp, N)                        # 4^al * wp = 1 mod 3^(N+1)
        for n in range(0, 600):
            v = v3int(2 ** n * w - 1)
            if (n - n0) % 2: assert v == 0; continue
            Mn = (n - n0) // 2
            pred = 1 + min(v3int(Mn - al) if Mn != al else 10 ** 9, N)
            if v <= N: assert v == pred, (w, n, v, pred)
            nchk += 1
        # best integer approximations of alpha_w below X = 3^20: max_{0<=M<3^20} v_3(M - alpha) - 20
        # = largest t with (alpha mod 3^t) < 3^20, minus 20
        # some M in [1, 3^20) with M = alpha mod 3^t  <=>  r := alpha mod 3^t satisfies 1 <= r < 3^20, or r = 0 and 3^t < 3^20
        tbest = max(t for t in range(0, N + 1) if (1 <= al % 3 ** t < 3 ** 20) or (al % 3 ** t == 0 and t < 20))
        rows.append((tbest - 20, w))
    print(f"  checked v_3(2^n w - 1) = [n=n0 mod 2](1 + v_3((n-n0)/2 - alpha_w)) for all w < 400 (3 !| w), n < 600:"
          f" {nchk} cases")
    rows.sort(reverse=True)
    print(f"  excess depth max_(1<=M<3^20) v_3(M - alpha_w) - 20 over w < 400: max {rows[0][0]} (w={rows[0][1]}),"
          f" distribution {sorted(set(r[0] for r in rows))}")
    print("  (w a power of 4 times 1 or 2: alpha_w is a non-positive integer and v_3 = 1 + v_3(M + a) exactly;"
          " otherwise alpha_w is a 3-adic irrational and the depth is the length of a digit coincidence)")
    for w in (1, 2, 4, 5, 7, 11, 13):
        n0 = 0 if w % 3 == 1 else 1
        al = dlog4((2 ** n0) * w, 30)
        digs = []; a = al
        for _ in range(20): digs.append(a % 3); a //= 3
        print(f"    alpha_{w:<3d} mod 3^20 digits (low first): {''.join(map(str, digs))}")

def section_L():
    print("=" * 100)
    print("L. Yu's p-adic bound (two logarithms, p = 3), explicit form")
    m = 2; p = 3
    C = (16 * math.e) ** (2 * (m + 1)) * m ** 1.5 * (math.log(2 * m)) ** 2 * p / (math.log(p)) ** 2
    print(f"  v_3(x1^b1 x2^b2 - 1) < C (log A1)(log A2) max(log T, delta B/B2),  C = (16e)^6 2^1.5 (log 4)^2 3/(log 3)^2"
          f" = {C:.4e}")
    print("  applied to 2^n w - p = p (2^(n-r) (2^r w/p) - 1), r in {0,1,2} with 3 !| n-r (so v_3 condition holds):")
    print("  v_3(2^n w - p) <= C * 1 * log(4 max(w,p,e)) * (log n + 7 log 3 + 17 + log log 4 max(w,p,e) + 2):"
          " O(log n log H)")
    for (n, w, p) in ((10 ** 3, 1, 1), (10 ** 6, 7, 1), (10 ** 9, 1, 43), (10 ** 18, 1000, 1)):
        H = max(w, p, math.e) * 4
        T = 2 * n / 0.5 * math.exp(3 * 17) * p ** 3 * 1.0     # 2B_m/delta e^{(m+1)(6m+5)} p^{m+1} log A1 (log A1 = 1)
        bound = C * 1.0 * math.log(H) * math.log(T)
        triv = n * math.log(2) / math.log(3) + math.log(w) / math.log(3)
        print(f"    n={n:.0e} w={w} p={p}: Yu bound {bound:.3e}   trivial bound n log_3 2 = {triv:.3e}"
              f"   ({'Yu better' if bound < triv else 'trivial better'})")

def section_M(AMAX=300, U=6000):
    print("=" * 100)
    print(f"M. double exhaustion: 2^(a+1) u - 1 = 3^k w' (1/2 thread) and 2^a u - 1 = 3^k w' (1 thread),"
          f" u odd <= {U}, w' <= {U}, k >= 3, a <= {AMAX}")
    sols = {'1/2': [], '1': []}
    for a in range(0, AMAX + 1):
        pa = 1 << a
        for u in range(1, U + 1, 2):
            if u % 3 == 0: continue
            for thr, t in (('1/2', 2 * pa * u - 1), ('1', pa * u - 1)):
                if t <= 0: continue
                k = v3int(t)
                if k >= 3:
                    wq = t // 3 ** k
                    if wq <= U: sols[thr].append((pa * u, a, u, k, wq))
    for thr in ('1/2', '1'):
        L = sorted(sols[thr])
        big = [s for s in L if s[0] > 10 ** 6]
        print(f"  {thr} thread: {len(L)} solutions; largest n = 2^a u: " +
              "; ".join(f"n=2^{s[1]}*{s[2]}={s[0]} = {'(1+3^%d*%d)/2' % (s[3], s[4]) if thr=='1/2' else '1+3^%d*%d' % (s[3], s[4])}"
                        for s in L[-4:]))
        print(f"    solutions with n > 10^6: {len(big)}; largest a among all: {max(s[1] for s in L)}")

def psi_int(m):
    t, half = (m - 1, False) if m % 3 == 1 else (2 * m - 1, True)
    k = 0
    while t % 3 == 0: t //= 3; k += 1
    return t << Kstar(k - 1), half, k

def section_N(WMAX=49, KMAX=400):
    print("=" * 100)
    print(f"N. digit-exhausted family m = (1 + 3^k w)/2, odd w <= {WMAX}, 3 !| w, 3 <= k <= {KMAX} (m up to ~10^{int(KMAX*0.477)+2})")
    worst = (0, None); worst_steps = (0, None); worstD = (0, None); cnt = 0
    lte_ok = 0; second = {}
    for w in range(1, WMAX + 1, 2):
        if w % 3 == 0: continue
        for k in range(3, KMAX + 1):
            m = (1 + 3 ** k * w) // 2
            x = m; mx = 1.0; steps = 0; DL1 = 0; Dt = 0
            while True:
                y, half, kk = psi_int(x); steps += 1; Dt += kk
                if half and kk >= 3: DL1 += kk
                if steps == 2: second.setdefault(w, []).append(kk)
                if w == 1 and steps == 1:
                    K = K0(k - 1); assert y == 1 << K
                    # LTE prediction for the next precision
                    if (K + 1) % 2 == 0: pred = 1 + v3int((K + 1) // 2)       # 2^(K+1) - 1 : 1/2 thread
                    else: pred = None
                    y2, h2, k2 = psi_int(y)
                    if h2: assert k2 == pred; lte_ok += 1
                    else: assert k2 == 1 + v3int(K // 2) and K % 2 == 0; lte_ok += 1
                x = y
                r = x / m
                if r > mx: mx = r
                if x < m: break
                assert steps < 100000
            cnt += 1
            if mx > worst[0]: worst = (mx, (w, k))
            if steps > worst_steps[0]: worst_steps = (steps, (w, k))
            if DL1 > worstD[0]: worstD = (DL1, (w, k))
    print(f"  {cnt} members: all descend under Psi; LTE predicted the second precision exactly in all {lte_ok} w=1 cases")
    print(f"  worst excursion {worst[0]:.4f} at (w,k)={worst[1]};  most Psi-steps {worst_steps[0]} at {worst_steps[1]};"
          f"  most L1 digits {worstD[0]} at {worstD[1]}")
    mx2 = max(max(v) for v in second.values())
    print(f"  precision of the second Psi-step (the landing 2^K0(k-1) w): max {mx2} over the family; "
          f"per w (max): " + ", ".join(f"{w}:{max(v)}" for w, v in sorted(second.items())[:12]))

def section_O():
    print("=" * 100)
    print("O. renewal structure of Psi under Haar measure on the 3-adic units")
    br = [(k, (Kstar(k - 1) + h) * LN2 - k * LN3) for h in (0, 1) for k in range(1, 400)]
    Ed = sum(k * 3.0 ** -k for k, l in br); El = sum(l * 3.0 ** -k for k, l in br)
    print(f"  branch (h,k) has Haar mass 3^-k, consumes k digits; E[digits] = {Ed:.6f}, E[ln rho] = {El:.6f},"
          f" drift {El/Ed:.6f} nats/digit (greedy G: ln(2/3) = {math.log(2/3):.6f})")
    Z = lambda s, th: sum(math.exp(-s * k * LN3 + th * l) for k, l in br)
    def s_of(th):
        lo, hi = 0.0, 1.5
        for _ in range(100):
            mid = (lo + hi) / 2
            if Z(mid, th) > 1: lo = mid
            else: hi = mid
        return lo
    lo_t, hi_t = 0.0, 6.0                                  # s(theta) is convex: ternary search
    for _ in range(80):
        m1 = lo_t + (hi_t - lo_t) / 3; m2 = hi_t - (hi_t - lo_t) / 3
        if s_of(m1) < s_of(m2): hi_t = m2
        else: lo_t = m1
    best = (s_of((lo_t + hi_t) / 2), (lo_t + hi_t) / 2)
    print(f"  dim_H(Bad_Psi) = min_theta s(theta), sum 3^(-s k) rho^theta = 1:  {best[0]:.5f} (theta = {best[1]:.3f})")
    print(f"  Haar-probability that Psi survives D digits ~ 3^(-{1-best[0]:.4f} D) = {3**-(1-best[0]):.4f}^D  (greedy G: 0.758751^D, dim 0.748)")
    f = lambda th: sum(3.0 ** -k * math.exp(th * l) for k, l in br)
    lo, hi = 1.0, 50.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 1: lo = mid
        else: hi = mid
    print(f"  Lundberg exponent theta* = {lo:.5f} (E[rho^theta*] = 1): random-model worst excursion over m <= X is X^{1/lo:.5f};"
          f"  X=1e11: {1e11**(1/lo):.2f}, X=1e18: {1e18**(1/lo):.2f}")
    @lru_cache(maxsize=None)
    def alive(n, A, D):
        tot = 2
        for half in (0, 1):
            for k in range(1, n + 1):
                A2 = A + Kstar(k - 1) + half; D2 = D + k
                if (1 << A2) < 3 ** D2: continue
                tot += alive(n - k, A2, D2)
        return tot
    line = []
    for n in (10, 20, 30, 37, 40, 60, 100, 150, 200):
        c = alive(n, 0, 0); line.append(f"{n}:{c} ({math.log(c, 3)/n:.4f})")
    print("  exact number of Psi-alive unit classes mod 3^(n+1) [n: count (log_3 count / n)]:")
    print("   " + "  ".join(line))
    print("  compare full choice (dimension lane, backward threads): 338 classes mod 3^42 (n=41), growth ~ n^1.7")

def main():
    sect = sys.argv[1] if len(sys.argv) > 1 else "KLMNO"
    random.seed(20260922)
    if "K" in sect: section_K()
    if "L" in sect: section_L()
    if "M" in sect: section_M()
    if "N" in sect: section_N()
    if "O" in sect: section_O()

if __name__ == '__main__':
    main()
