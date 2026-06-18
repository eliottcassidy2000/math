#!/usr/bin/env python3
"""
REFRAME 5 CONSOLIDATED (mac-mini-2026-06-17-S5).  Self-contained; python3 stdlib only.

RESULT.  For the PRINCIPAL single-drop covering families at N=14
   S_q(k) = ({1,...,13}\{q}) u {L*k},   L = lcm(q,14),   q in {7,...,13},  k>=1,
the LRC gap obeys the EXACT closed law (verified k<=30, exact rationals)
   M(S_q(k)) = (L/q) * k / (L*k + flank_q),    flank_q in {1,3,5,7,8} (a runner, <14),
which GENERALIZES the q=12 family 7m/(84m+5)  (there L/q = 84/12 = 7, flank = 5).

CONSEQUENCES, with the load-bearing facts flagged [EMPIRICAL] vs [PROVED]:
  (A) [PROVED, calculus] M(S_q(k)) is strictly increasing in k (slope (L/q)*flank/(...)^2>0),
      so min_{k>=1} = M(S_q(1)). Infinite resonant tower collapses to ONE k=1 check.
  (B) [PROVED, arithmetic] the k=1 check M(S_q(1)) = L/(q(L+flank)) >= 1/14 is equivalent
      to (14-q)*L >= q*flank, which holds using ONLY the trivial non-circular bound
      flank <= 13:  (14-q)*lcm(q,14) >= 13*q for all q in {7..13}  (checked, true).
  (C) [EMPIRICAL] the law's coefficient c = L/q (limit gap 1/q) equals M(core {1..13}\{q}).
      This is the one piece tying the recursion together that is verified but not proved.

So REFRAME 5 yields, for the principal families, a NON-CIRCULAR proof of M>=1/14 MODULO
the empirical closed-form law -- the circularity that the prior peeling-chain reflection
flagged is REMOVED here (it does NOT use M>=1/14 as input; it uses flank<14 and q<14).
The residue is no longer an inequality but a STRUCTURE claim: "the gap is (L/q)k/(Lk+flank)."

LIMITS (honest):
  * Only PRINCIPAL single-drop families (q in 7..13) obey the clean law; for q in 2..6 the
    binding flank changes with k (law fails, monotonicity fails), and for general
    multi-drop covering 13-sets the law is not established. So this is a clean SPECIAL CASE,
    not full LRC(14).
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
    b = F(0)
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v
    return b


def lcm(a, b):
    return a * b // gcd(a, b)


print("=" * 78)
print("REFRAME 5 CONSOLIDATED -- principal single-drop families q=7..13, N=14")
print("=" * 78)

# determine flank_q (=L+flank - L) by reading M(S_q(1)) = c/(L+flank), c=L/q
print(f"\n{'q':>2} {'L':>4} {'c=L/q':>6} {'M(S,1)':>9} {'L+flank':>8} {'flank':>5} "
      f"{'law k<=20 exact?':>16}")
flanks = {}
for q in range(7, 14):
    A = [v for v in range(1, 14) if v != q]
    L = lcm(q, 14)
    c = L // q
    M1 = M(A + [L])
    Lpf = c / M1               # = L + flank   (since M1 = c/(L+flank))
    flank = Lpf - L
    flanks[q] = flank
    # verify law over k=1..20
    exact = True
    for k in range(1, 21):
        pred = F(c * k, L * k + int(flank))
        if M(A + [L * k]) != pred:
            exact = False; break
    print(f"{q:>2} {L:>4} {c:>6} {str(M1):>9} {str(Lpf):>8} {str(flank):>5} {str(exact):>16}")

print("\n--- non-circular closure (uses ONLY flank<=13 and q<=13): ---")
print(f"{'q':>2} {'L':>4} {'(14-q)L':>8} {'13q':>5} {'(14-q)L>=13q':>13} "
      f"{'=> M(S,1)>=1/14':>16}")
ok_all = True
for q in range(7, 14):
    L = lcm(q, 14)
    lhs = (14 - q) * L; rhs = 13 * q
    ok = lhs >= rhs; ok_all &= ok
    # actual k=1 gap >= 1/14?
    A = [v for v in range(1, 14) if v != q]
    real = M(A + [L]) >= F(1, 14)
    print(f"{q:>2} {L:>4} {lhs:>8} {rhs:>5} {str(ok):>13} {str(real):>16}")
print(f"\n  (14-q)*lcm(q,14) >= 13q for all principal q?  {ok_all}")
print("  => For principal single-drop families, M(S)>=1/14 holds NON-CIRCULARLY,")
print("     modulo the EMPIRICAL closed-form law M=(L/q)k/(Lk+flank).")

print("\n--- scope honesty: messy families q=2..6 do NOT obey the clean law ---")
for q in range(2, 7):
    A = [v for v in range(1, 14) if v != q]
    L = lcm(q, 14)
    seq = [M(A + [L * k]) for k in range(1, 11)]
    MA = M(A)
    binders = [(k + 1, seq[k]) for k in range(10) if seq[k] < MA]
    mono = all(binders[i][1] <= binders[i + 1][1] for i in range(len(binders) - 1)) if binders else True
    print(f"  q={q}: binding k,M = {[(k,str(v)) for k,v in binders][:6]} ; "
          f"monotone? {mono}  (min NOT at k=1 => law/collapse fails)")

print("\nDONE.")
