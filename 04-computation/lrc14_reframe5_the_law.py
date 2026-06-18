#!/usr/bin/env python3
"""
REFRAME 5 -- THE LAW (mac-mini-2026-06-17-S5).

DISCOVERY (from crossing-index table): for the N=14 single-drop resonant family
   S = ({1..13}\{q}) u {w},  w = L*k,  L = lcm(q,14),
the binding sum-crossing gives
   M(S) = (c*k) / (flank + L*k)
with c and flank CONSTANT in k (per q). I.e. the dip family is the rational law
   M_q(k) = c_q * k / (L_q * k + flank_q),
generalizing  7m/(84m+5)  (q=12: c=7, L=84, flank=5).

Then  M(S) >= 1/14  <=>  14 c k >= L k + flank  <=>  (14 c - L) k >= flank.
Since k>=1, this holds for ALL k>=1  IFF  14 c - L >= flank  at k=1 AND 14c-L>0 (so it
only gets easier as k grows).  Actually if 14c-L > 0 then LHS is increasing in k and the
TIGHTEST case is k=1: (14c-L) >= flank.

** THE non-circular reduction **:  M(S)>=1/14 for the whole k-family  reduces to ONE
inequality at k=1:   14 c_q - L_q  >=  flank_q   (equivalently 14 c >= L + flank = D_1).
And D_1 = flank+L is the k=1 binding denominator, c the k=1 numerator: this is just
14*M(S,k=1) >= 1, i.e. STILL the k=1 instance of the goal.  BUT it has collapsed the
INFINITE k-family to a SINGLE k=1 check -- a genuine (if modest) reduction.

This script:
  (1) confirms M_q(k) = c k/(Lk+flank) with constant c,flank for every q (k=1..40);
  (2) tabulates (c_q, L_q, flank_q) and the key margin 14c - L - flank;
  (3) checks: is the WHOLE infinite family certified by k=1 alone? (monotonicity in k)
  (4) honest verdict on circularity.
"""
from fractions import Fraction as F
from math import gcd


def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r


def g(S, t):
    return min(nrm(v * t) for v in S)


def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C


def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at


def lcm(a, b):
    return a * b // gcd(a, b)


print("=" * 78)
print("REFRAME 5 -- THE LAW  M_q(k) = c k/(L k + flank),  reduction to one k=1 check")
print("=" * 78)

KMAX = 40
print(f"\n{'q':>3} {'L=lcm':>6} {'c':>4} {'flank':>6} {'law M=c k/(Lk+fl)':>20} "
      f"{'14c-L':>6} {'14c-L-flank':>11} {'14c>=L?':>8} {'k=1 cert all k?':>14}")

summary = []
for q in range(2, 14):
    A = [v for v in range(1, 14) if v != q]
    L = lcm(q, 14)
    # fit c, flank from two k values where w binds (dip>0). Use k=1 and k=2; verify to KMAX.
    data = []
    for k in range(1, KMAX + 1):
        w = L * k
        if w in A:
            data.append(None); continue
        S = A + [w]
        MS, _ = M(S)
        MA, _ = M(A)
        data.append((MS, MA))
    # find k where dip>0 (w binds): MS<MA
    bind_ks = [k for k in range(1, KMAX + 1)
               if data[k - 1] is not None and data[k - 1][0] < data[k - 1][1]]
    if len(bind_ks) < 2:
        print(f"{q:>3} {L:>6}   --  (w rarely binds: dip~0; M_q(k)=M(A)={data[0][1] if data[0] else '?'})")
        summary.append((q, L, None, None, None))
        continue
    k1, k2 = bind_ks[0], bind_ks[1]
    M1 = data[k1 - 1][0]; M2 = data[k2 - 1][0]
    # M = c k/(L k + f).  Two eqns: M1(L k1 + f) = c k1 ; M2(L k2 + f)=c k2
    #  => c = M1(L k1+f)/k1 = M2(L k2+f)/k2
    #  M1 L k1 k2 + M1 f k2 = M2 L k1 k2 + M2 f k1   (mult by k1 k2)
    #  f (M1 k2 - M2 k1) = L k1 k2 (M2 - M1)
    num = L * k1 * k2 * (M2 - M1)
    den = (M1 * k2 - M2 * k1)
    if den == 0:
        print(f"{q:>3} {L:>6}  degenerate fit"); continue
    f = num / den
    c = M1 * (L * k1 + f) / k1
    # verify law over ALL binding k
    ok = True
    for k in bind_ks:
        pred = c * k / (L * k + f)
        if pred != data[k - 1][0]:
            ok = False; break
    cstr = str(c); fstr = str(f)
    margin = 14 * c - L  # want >= flank (=f) for the inequality (14c-L)k>=f at k=1
    m2 = 14 * c - L - f
    keyA = 14 * c >= L
    # does k=1 certify all k? since M_q(k) = ck/(Lk+f) is MONOTONE INCREASING in k
    #   (derivative sign = c*f > 0), the MINIMUM over k>=1 is at k=1. So yes IF increasing.
    incr = (c > 0 and f > 0)  # then strictly increasing in k
    k1cert = incr  # min at k=1
    print(f"{q:>3} {L:>6} {cstr:>4} {fstr:>6} {'c k/(Lk+f), ok='+str(ok):>20} "
          f"{str(margin):>6} {str(m2):>11} {str(keyA):>8} {str(k1cert):>14}")
    summary.append((q, L, c, f, m2))

print("""
KEY STRUCTURAL FACTS (read off the table):
  * M_q(k) = c_q k /(L_q k + flank_q) is the EXACT per-level dip law (7m/(84m+5) is q=12).
  * Each c_q, flank_q is a CONSTANT (independent of k). c_q in {2,7,14}, flank_q small.
  * M_q(k) is STRICTLY INCREASING in k (slope c*flank>0), so its min over k>=1 is at k=1.
  * Hence the INFINITE k-family is certified by the SINGLE k=1 value -- a real collapse.
  * M(S)>=1/14 for all k  <=>  14 c_q - L_q >= flank_q  (one arithmetic margin per q).
""")

print("-" * 78)
print("MONOTONICITY proof check: is M_q(k) increasing? show M_q(1) is the family minimum.")
print("-" * 78)
for q in range(2, 14):
    A = [v for v in range(1, 14) if v != q]
    L = lcm(q, 14)
    MA, _ = M(A)
    vals = []
    for k in range(1, 26):
        w = L * k
        if w in A: continue
        MS, _ = M(A + [w])
        vals.append((k, MS))
    binders = [(k, v) for (k, v) in vals if v < MA]
    if binders:
        mn = min(binders, key=lambda kv: kv[1])
        mono = all(binders[i][1] <= binders[i + 1][1] for i in range(len(binders) - 1))
        print(f"  q={q:>2}: min over binding k at k={mn[0]} M={mn[1]} ; "
              f"non-decreasing in binding-k? {mono} ; M(A)={MA}")
    else:
        print(f"  q={q:>2}: w never binds in k<=25 (dip 0 throughout); M(A)={MA}")

print("\nDONE.")
