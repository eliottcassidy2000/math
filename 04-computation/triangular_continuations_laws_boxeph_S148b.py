#!/usr/bin/env python3
"""
triangular_continuations_laws_boxeph_S148b.py  (HYP-8175 verification layer)

The wide-birth scan (S148a) surfaced candidate laws; this script PROVES/verifies
each exactly and factors the recurrences.

 LAW A (SHEAR-DELAY LAW, proved): for the subsets grid C(n,x), the shear-t row
   sums A_t(m) = sum_x C(m - t*x, x) satisfy  A_t(m) = A_t(m-1) + A_t(m-t-1).
   Proof: Pascal C(m-tx, x) = C(m-1-tx, x) + C((m-t-1)-t(x-1), x-1); resum.  QED.
   t=0: 2^m; t=1: Fibonacci; t=2: Narayana's cows; t=3: A003269 (delay-4), ...
   THE OWNER'S "moved down t" = delay-(t+1) generalized Fibonacci.  Verified to t=6.

 LAW B (MOSER-SHEAR = FIBONACCI MINUS ONE, proved): for the truncation grid
   R_x(n) = sum_{i<=x} C(n,i), the shear-1 sums equal F(m+3) - 1.
   Proof: swap sums, hockey-stick sum_{n<=m-i} C(n,i) = C(m-i+1, i+1), then the
   shallow-diagonal identity sum_j C(m+2-j, j) = F(m+3); the j=0 term is 1.  QED.

 LAW C (q-CONFLUENCE, proved): [n,2]_q = (q^{2(n-1)+...})... has char poly
   (x-1)(x-q)(x-q^2): the q-triangular numbers are the {1,q,q^2}-geometric
   mixture; q->1 is the confluent triple root => polynomial degree 2 (ordinary
   triangular).  Explicit partial fractions verified for q=2..7.

 LAW D (OWNER ROW SUMS = FAULHABER SHEAR-1 + CASCADE): rowsum_T(m) =
   [sum_x S_x(m-x)] + 2^{m-6} for m >= 6 (and equal below).  Exact restatement
   of the S147 closed form summed over j; verified m <= 16 (rows >= 8 assume the
   S147 form, which is verified only on the given 28 entries -- honest caveat).
   Similarly ownerT1(m) = faulhaberT2(m) + cascade (first split at m = 9).

 LAW E (PROTH SHEAR-2 CHARACTERISTIC): the r=5 recurrence found by the detector
   factors as (x-1)^2 (x+1) (x^2-2): the sqrt(2)-geometric pair (from 2^{m/2})
   + parity + the linear-in-K drift.  Matches the hand closed form.
   Also moser shear-2's r=5 [1,1,0,0,-1] factors as (x^2-x-1)(x^3-...)? -- computed.

boxeph-2026-07-20-S148.  Exact.
"""

from fractions import Fraction as Fr
from math import comb

# ---------------------------------------------------------------- LAW A
print("=" * 96)
print("LAW A: shear-delay law  A_t(m) = A_t(m-1) + A_t(m-t-1)   (subsets grid)")
print("=" * 96)
def A_t(t, m): return sum(comb(m - t * x, x) for x in range(0, m + 1) if m - t * x >= x >= 0)
names = {0: "2^m", 1: "Fibonacci", 2: "Narayana's cows A000930", 3: "A003269", 4: "A003520",
         5: "A005708", 6: "A005709(delay-7)"}
for t in range(0, 7):
    seq = [A_t(t, m) for m in range(0, 18)]
    ok = all(seq[m] == seq[m - 1] + seq[m - t - 1] for m in range(t + 1, 18))
    print("  t=%d: %s...  recurrence a(m)=a(m-1)+a(m-%d): %s   [%s]" %
          (t, seq[:9], t + 1, "OK" if ok else "FAIL", names.get(t, "")))
    assert ok
print("  PROOF: Pascal on the top index; the two terms re-index to (m-1) and (m-t-1). QED")

# ---------------------------------------------------------------- LAW B
print("\n" + "=" * 96)
print("LAW B: moser(trunc) shear-1 sums = F(m+3) - 1")
print("=" * 96)
def R(x, n): return sum(comb(n, i) for i in range(0, x + 1))
fib = [0, 1]
while len(fib) < 40: fib.append(fib[-1] + fib[-2])
mos = [sum(R(x, m - x) for x in range(0, m + 1)) for m in range(0, 24)]
ok = all(mos[m] == fib[m + 3] - 1 for m in range(0, 24))
print("  sums:", mos[:10], " = F(m+3)-1:", ok)
assert ok
print("  PROOF: swap sums; hockey-stick sum_{n<=m-i}C(n,i) = C(m-i+1,i+1); shallow")
print("  diagonals sum to F(m+3); subtract the j=0 term.  QED")

# ---------------------------------------------------------------- LAW C
print("\n" + "=" * 96)
print("LAW C: q-confluence -- [n,2]_q has char poly (x-1)(x-q)(x-q^2); q=1 = triple root")
print("=" * 96)
def qb2(q, n):
    if n < 2: return 0
    if q == 1: return comb(n, 2)
    return (q ** n - 1) * (q ** (n - 1) - 1) // ((q * q - 1) * (q - 1))
for q in range(2, 8):
    seq = [qb2(q, n) for n in range(2, 14)]
    c1, c2, c3 = 1 + q + q * q, -(q + q * q + q ** 3), q ** 3
    ok = all(seq[i] == c1 * seq[i - 1] + c2 * seq[i - 2] + c3 * seq[i - 3] for i in range(3, len(seq)))
    # partial fractions: [n,2]_q = (q^{2n} ... ) explicit mixture a*1 + b*q^n + c*q^{2n}
    D = (q * q - 1) * (q - 1) * q
    aa = Fr(1, (q - 1) * (q * q - 1))          # coefficient of 1 (n-free)
    bb = Fr(-1, (q - 1) * (q - 1) * q)         # coefficient of q^n ((q+1) cancels vs q^2-1)
    cc = Fr(1, (q * q - 1) * (q - 1) * q)      # coefficient of q^{2n}
    mix = all(Fr(qb2(q, n)) == aa + bb * q ** n + cc * q ** (2 * n) for n in range(2, 12))
    print("  q=%d: rec (1+q+q^2, -q(1+q+q^2), q^3) OK=%s ; mixture a+b q^n+c q^{2n} with" % (q, ok))
    print("        a=%s b=%s c=%s : %s" % (aa, bb, cc, mix))
    assert ok and mix
print("  => triangular = confluent limit of three geometric strands {1, q, q^2}.")

# ---------------------------------------------------------------- LAW D
print("\n" + "=" * 96)
print("LAW D: owner row sums = faulhaber-shear1 + 2^(m-6) cascade  (m>=6)")
print("=" * 96)
def S(p, n): return sum(k ** p for k in range(1, n + 1))
def T_owner(m, j):
    c = comb(m - 6, j - 4) if 0 <= j - 4 <= m - 6 else 0
    return S(j - 1, m - j + 1) + c
rows = {m: [T_owner(m, j) for j in range(1, m + 1)] for m in range(1, 17)}
given = {1:[1],2:[2,1],3:[3,3,1],4:[4,6,5,1],5:[5,10,14,9,1],6:[6,15,30,37,17,1],7:[7,21,55,101,99,33,1]}
assert all(rows[m] == given[m] for m in given)
fT1 = [sum(S(x, m - x) for x in range(0, m + 1)) for m in range(0, 17)]
ok6 = all(sum(rows[m]) == fT1[m] + (2 ** (m - 6) if m >= 6 else 0) for m in range(1, 17))
print("  rowsum_T(m) = fT1(m) + [m>=6] 2^(m-6):", ok6,
      " e.g. m=6: %d = %d + 1 ; m=7: %d = %d + 2" % (sum(rows[6]), fT1[6], sum(rows[7]), fT1[7]))
assert ok6
fT2 = [sum(S(x, m - 2 * x) for x in range(0, m // 2 + 1)) for m in range(0, 17)]
oT1 = [sum(T_owner(m - x, x + 1) for x in range(0, m) if 1 <= x + 1 <= m - x) for m in range(1, 17)]
split = [(m, oT1[m - 1] - fT2[m]) for m in range(1, 17)]
print("  ownerT1(m) - faulhaberT2(m):", split[:12], "(cascade joins at m=9: C(m-x-6,x-3))")
print("  CAVEAT: rows >= 8 of T are the S147 closed form (verified only on the 28 given")
print("  entries); every law above is exact for the closed-form extension.")

# ---------------------------------------------------------------- LAW E
print("\n" + "=" * 96)
print("LAW E: characteristic factorizations of the detected recurrences")
print("=" * 96)
def polydiv(num, den):
    num = num[:]
    out = [0] * (len(num) - len(den) + 1)
    while len(num) >= len(den):
        f = Fr(num[-1], den[-1]); d = len(num) - len(den)
        out[d] = f
        for i in range(len(den)): num[d + i] -= f * den[i]
        while num and num[-1] == 0: num.pop()
    return out, num
# proth t=2: r=5 coeffs [1,3,-3,-2,2] => char x^5 - x^4 - 3x^3 + 3x^2 + 2x - 2
P = [-(-2), -2, 3, -3, -1, 1]     # ascending: constant..x^5  of x^5 -x^4 -3x^3 +3x^2 +2x -2
P = [-2, 2, 3, -3, -1, 1]
for root, tag in ((1, "(x-1)"), (1, "(x-1) again"), (-1, "(x+1)")):
    q, r = polydiv(P, [-root, 1])
    if not r: P = [Fr(v) for v in q]; print("  proth-t2 char / %s exact" % tag)
print("  remaining factor coeffs (ascending):", P, " = x^2 - 2 :", P == [Fr(-2), Fr(0), Fr(1)])
# moser t=2: [1,1,0,0,-1] => x^5 - x^4 - x^3 + 1
Q = [1, 0, 0, -1, -1, 1]
q, r = polydiv(Q, [-1, 1])
print("  moser-t2 char / (x-1):", "exact" if not r else "NOT a factor", "->", q if not r else "")
if not r:
    q2, r2 = polydiv(q, [-1, -1, 1])   # x^2 - x - 1
    print("    then / (x^2-x-1):", "exact -> Fibonacci strand present! rest=%s" % q2 if not r2
          else "not a factor (rest %s)" % r2)
print("\nDONE.")
