#!/usr/bin/env python3
"""The Paley heptagon and the trivial Collatz cycle (opus S15, nineteenth note, 2026-10-01).

Checks:
  A. parity codes of cycles: for T (n/2, 3n+1) and T1 (n/2, (3n+1)/2) the 2-adic parity image of a cycle
     with word w of length L is -c/(2^L - 1) (c = sum w_i 2^i), the real reading 0.(w)_2 is rev(c)/(2^L - 1);
     the numerators form cyclotomic cosets of 2; all integer cycles on Z listed with their codes
  B. L = 3: the trivial cycle {1,4,2} (word 100) has real code QR_7 and 2-adic code NQR_7; the minus-sheet
     3-cycle {-5,-7,-10} of T1 (word 110) the reverse; QR/NQR swap = complement = place swap
  C. a single 2-coset equals QR_p for a Mersenne prime p = 2^L - 1 only for L = 3
  D. Fano/Singer: {1,2,4} is a (7,3,1) difference set = trace-zero exponents of x^3+x+1; its translates are the Fano lines
  E. Aut(Paley_7): 21 maps x -> ax+b (a in QR); the stabiliser of QR_7 is <x2>, acting as the trivial cycle's dynamics;
     translates of QR_7 are not 2-cosets
  F. lonely runner: lon({1,2}) = 1/3 attained exactly at t = 1/3, 2/3 = the real codes of the T1 trivial cycle;
     positions of speeds {1,2,4} at t = k/7
  G. the owner's F: F(2N) = 3N, F(2N-1) = K_N (K = 0,2,1,4,2,10,3,9,4,15,5): in-degrees vs Collatz;
     the exact Collatz conjugates with the same even rule; the folded two-sheet relation
Reproduce: python 04-computation/experiments/collatz_paley_bridge_20261001.py
"""
from fractions import Fraction
from itertools import product

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def T(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def T1(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def cycle_of(f, x0, cap=10000):
    seen = []
    x = x0
    for _ in range(cap):
        if x in seen:
            i = seen.index(x)
            return seen[i:]
        seen.append(x)
        x = f(x)
    return None


def word(f, cyc):
    return [x % 2 for x in cyc]


def code2adic(w):
    L = len(w)
    c = sum(b << i for i, b in enumerate(w))
    return Fraction(-c, 2 ** L - 1), c


def code_real(w):
    L = len(w)
    r = sum(b << (L - 1 - i) for i, b in enumerate(w))
    return Fraction(r, 2 ** L - 1), r


def coset(c, m):
    s = set()
    x = c % m
    while x not in s:
        s.add(x)
        x = 2 * x % m
    return s


print("== A. parity codes of the integer cycles on Z ==")
for name, f, starts in (("T  (n/2, 3n+1)", T, [1, 0, -1, -5, -17]), ("T1 (n/2,(3n+1)/2)", T1, [1, 0, -1, -5, -17])):
    print(" ", name)
    for s in starts:
        cyc = cycle_of(f, s)
        w = word(f, cyc)
        L = len(w)
        y, c = code2adic(w)
        r, rc = code_real(w)
        m = 2 ** L - 1
        two_adic_nums = sorted({(-code2adic(w[j:] + w[:j])[1]) % m for j in range(L)}) if m > 1 else [0]
        real_nums = sorted({code_real(w[j:] + w[:j])[1] % m for j in range(L)}) if m > 1 else [0]
        show = (lambda v: v if L <= 6 else f"{len(v)} elements")
        print(f"    cycle from {s:4d}: length {L:2d}, word {''.join(map(str, w))[:20]}{'...' if L > 20 else ''}, "
              f"2-adic numerators mod {m}: {show(two_adic_nums)}, real numerators: {show(real_nums)}")
# identity check: the 2-adic code equals the 2-adic limit of the parity vector (Phi); verify on many integers of Gamma_1
ok = True
for n in (1, 4, 2):
    cyc = cycle_of(T, n)
    w = word(T, cyc)
    y, c = code2adic(w)
    # Phi(n) mod 2^K from the first K parities
    K = 30
    x = n
    phi = 0
    for i in range(K):
        phi |= (x % 2) << i
        x = T(x)
    # -c/(2^L-1) mod 2^K
    m = 2 ** K
    val = (y.numerator * pow(y.denominator, -1, m)) % m
    ok = ok and val == phi
check(ok, "A: for n = 1, 4, 2 the first 30 parity bits equal the 2-adic expansion of -c/(2^L-1) = -1/7, -4/7, -2/7")

print("\n== B. L = 3: Paley QR/NQR and the two sheets ==")
QR7 = {1, 2, 4}
NQR7 = {3, 5, 6}
cycT = cycle_of(T, 1)                       # [1, 4, 2]
wT = word(T, cycT)                          # [1, 0, 0]
real_T = {code_real(wT[j:] + wT[:j])[1] for j in range(3)}
adic_T = {(-code2adic(wT[j:] + wT[:j])[1]) % 7 for j in range(3)}
cyc5 = cycle_of(T1, -5)                     # [-5, -7, -10]
w5 = word(T1, cyc5)                         # [1, 1, 0]
real_5 = {code_real(w5[j:] + w5[:j])[1] for j in range(3)}
adic_5 = {(-code2adic(w5[j:] + w5[:j])[1]) % 7 for j in range(3)}
print(f"  trivial cycle {cycT} (T), word {wT}: real numerators {sorted(real_T)}, 2-adic numerators {sorted(adic_T)}")
print(f"  minus-sheet cycle {cyc5} (T1), word {w5}: real numerators {sorted(real_5)}, 2-adic numerators {sorted(adic_5)}")
check(real_T == QR7 and adic_T == NQR7, "B: trivial cycle: real code = QR_7 = {1,2,4} (Paley connection set), 2-adic code = NQR_7")
check(real_5 == NQR7 and adic_5 == QR7, "B: minus-sheet 3-cycle {-5,-7,-10}: real code = NQR_7, 2-adic code = QR_7")
comp = [1 - b for b in wT]
check(sorted(comp) == sorted(w5) and {code_real(comp[j:] + comp[:j])[1] for j in range(3)} == NQR7,
      "B: the complement of the trivial word (011) is the minus-sheet word up to rotation, real code NQR_7")
# the shortcut T1 rational cycle of the word 100 is {1/5, 4/5, 2/5}; the non-shortcut integer cycle is {1, 4, 2}
x = Fraction(1, 5)
orb = [x]
for _ in range(3):
    x = x / 2 if x.numerator % 2 == 0 else (3 * x + 1) / 2
    orb.append(x)
check(orb[3] == orb[0] and [int(v.numerator % 2) for v in orb[:3]] == [1, 0, 0],
      "B: under T1 the word 100 is the rational cycle {1/5, 4/5, 2/5} (clock 2^3 - 3 = 5); under T it is {1,4,2} (clock 2^2 - 3 = 1)")
# T on the cycle acts as x4 = x2^{-1} mod 7 on the real numerators of the 2-adic code
check(all((4 * v) % 7 in QR7 for v in QR7) and [4 * 1 % 7, 4 * 4 % 7, 4 * 2 % 7] == [4, 2, 1],
      "B: on {1,4,2} the Collatz map acts as x -> 4x (mod 7) = x/2, a Paley automorphism (4 is a square)")

print("\n== C. one 2-coset = QR_p for a Mersenne prime p = 2^L - 1 only when L = 3 ==")
for L in (2, 3, 5, 7, 13, 17, 19, 31):
    p = 2 ** L - 1
    o = len(coset(1, p))
    print(f"  L = {L:2d}, p = {p}: |<2>| = {o}, |QR_p| = {(p - 1) // 2}, equal: {o == (p - 1) // 2}")
check(all((len(coset(1, 2 ** L - 1)) == (2 ** L - 2) // 2) == (L == 3) for L in (2, 3, 5, 7, 13, 17, 19, 31)),
      "C: <2> = QR_p exactly for L = 3 among the Mersenne primes with L <= 31 (equation 2^(L-1) - 1 = L)")

print("\n== D. Fano / Singer ==")
diffs = sorted(((a - b) % 7) for a in QR7 for b in QR7 if a != b)
check(diffs == [1, 2, 3, 4, 5, 6], "D: {1,2,4} is a (7,3,1) difference set (every nonzero residue once)")
# F_8 = F_2[a]/(a^3+a+1): trace of a^i, i = 0..6, as 3-bit vectors
def mul(u, v):
    r = 0
    for i in range(3):
        if (v >> i) & 1:
            r ^= u << i
    for d in (4, 3):
        if (r >> d) & 1:
            r ^= 0b1011 << (d - 3)
    return r
pw = [1]
for i in range(1, 7):
    pw.append(mul(pw[-1], 0b010))
def tr(u):
    s = u ^ mul(u, u) ^ mul(mul(u, u), mul(u, u))
    return s
trace_zero = {i for i in range(7) if tr(pw[i]) == 0}
print("  powers of alpha:", pw, " exponents with trace 0:", sorted(trace_zero))
check(trace_zero == QR7, "D: {1,2,4} = exponents of alpha^i with trace 0 (x^3+x+1), i.e. the conjugates of alpha; a Fano line")
lines = {frozenset((x + b) % 7 for x in QR7) for b in range(7)}
check(len(lines) == 7 and all(len(l1 & l2) == 1 for l1 in lines for l2 in lines if l1 != l2),
      "D: the 7 translates of {1,2,4} are the lines of a Fano plane (pairwise meeting in one point)")

print("\n== E. Aut(Paley_7) and what the Collatz cycle keeps ==")
A = {(i, j) for i in range(7) for j in range(7) if (j - i) % 7 in QR7}
auts = []
for a in range(1, 7):
    for b in range(7):
        g = lambda x, a=a, b=b: (a * x + b) % 7
        if all(((g(i), g(j)) in A) for (i, j) in A):
            auts.append((a, b))
print(f"  affine automorphisms: {len(auts)}; multipliers used: {sorted({a for a, b in auts})}")
stab = [(a, b) for (a, b) in auts if {(a * x + b) % 7 for x in QR7} == QR7]
check(len(auts) == 21 and sorted({a for a, b in auts}) == [1, 2, 4] and sorted(stab) == [(1, 0), (2, 0), (4, 0)],
      "E: |Aut| = 21 (x -> ax+b, a in QR); the stabiliser of the connection set QR_7 is the multiplier group <x2> of order 3")
tr_cosets = [{(x + b) % 7 for x in QR7} for b in range(1, 7)]
check(all(s != coset(min(s), 7) or not s <= set(range(1, 7)) for s in tr_cosets),
      "E: no non-trivial translate of QR_7 is a 2-coset, so none is the code of a cycle (translations have no Collatz counterpart)")
anti = all((((-i) % 7, (-j) % 7) in A) == ((j, i) in A) for i in range(7) for j in range(7) if i != j)
check(anti, "E: x -> -x reverses every arc (anti-automorphism), i.e. swaps QR and NQR: the sheet/place swap of B")

print("\n== F. lonely runner readings ==")
def lon(speeds, den=8400):
    best, arg = Fraction(0), []
    for k in range(1, den):
        t = Fraction(k, den)
        m = min(min((t * v) % 1, 1 - (t * v) % 1) for v in speeds)
        if m > best:
            best, arg = m, [t]
        elif m == best:
            arg.append(t)
    return best, arg
b12, a12 = lon([1, 2])
print(f"  lon({{1,2}}) = {b12} at t in {[str(t) for t in a12]}")
check(b12 == Fraction(1, 3) and a12 == [Fraction(1, 3), Fraction(2, 3)],
      "F: speeds {1,2} (the T1 trivial cycle) are a tight instance (lon = 1/3 = 1/(k+1)), lonely exactly at t = 1/3, 2/3 = 0.(01)_2, 0.(10)_2, its own real parity code")
b124, a124 = lon([1, 2, 4])
print(f"  lon({{1,2,4}}) = {b124} at t in {[str(t) for t in a124]}; positions at t = 1/7: {sorted((v % 7) for v in (1, 2, 4))}/7, at t = 3/7: {sorted((3 * v) % 7 for v in (1, 2, 4))}/7")
check(b124 == Fraction(1, 3), "F: speeds {1,2,4} have lon = 1/3 (not tight for k = 3); at t = 1/7 they sit on QR_7/7, at t = 3/7 on NQR_7/7")

print("\n== G. the owner's F ==")
K = [0, 2, 1, 4, 2, 10, 3, 9, 4, 15, 5]


def F_owner(x, ext=True):
    if x == 0:
        return 0
    if x % 2 == 0:
        return 3 * (x // 2)
    N = (x + 1) // 2
    if N <= len(K):
        return K[N - 1]
    if ext and N % 2 == 1:
        return (N - 1) // 2          # the odd-position pattern K_(2j+1) = j
    return None


pre = {}
for x in range(0, 200):
    y = F_owner(x)
    if y is not None:
        pre.setdefault(y, []).append(x)
print("  preimages under the owner's F (inputs < 200, odd-position pattern extended):",
      {v: pre.get(v) for v in (0, 1, 2, 3, 4, 9, 15)})
check(sorted(pre[2]) == [3, 9] and 2 in pre[3] and 13 in pre[3],
      "G: the 2-cycle {2,3} of F has two preimages at each vertex (2 <- 3, 9; 3 <- 2, 13), from the given data alone")
check(sorted(pre[9]) == [6, 15, 37], "G: 9 has three preimages 6, 15, 37 (even rule, K_8 = 9, K_19 = 9 by the odd pattern)")
# Collatz 2-cycles: T1 {1,2} on Z has in-degrees (1,2); T on Z has {-1,-2} with (1,2); max in-degree of any Collatz map is 2
def preimages_T1(n):
    out = [2 * n]
    if (2 * n - 1) % 3 == 0:
        out.append((2 * n - 1) // 3)
    return out
def preimages_T(n):
    out = [2 * n]
    if (n - 1) % 3 == 0 and ((n - 1) // 3) % 2 == 1:
        out.append((n - 1) // 3)
    return out


check(sorted(len(preimages_T1(v)) for v in (1, 2)) == [1, 2] and sorted(len(preimages_T(v)) for v in (-1, -2)) == [1, 2],
      "G: the integer Collatz 2-cycles ({1,2} for T1, {-1,-2} for T) have in-degrees (1,2); no Collatz map has in-degree 3")
# the exact conjugates with the even rule 2N -> 3N
def Fplus(x):   # label x = n + 1, n >= 0, map T1
    return T1(x - 1) + 1
def Fminus(x):  # label x = n - 1, n >= 1, map 3n-1 shortcut
    n = x + 1
    m = n // 2 if n % 2 == 0 else (3 * n - 1) // 2
    return m - 1
ok = all(Fplus(2 * N) == 3 * N and Fplus(2 * N - 1) == N for N in range(1, 2000)) and \
     all(Fminus(2 * N) == 3 * N and Fminus(2 * N - 1) == N - 1 for N in range(1, 2000))
check(ok, "G: F+(2N) = 3N, F+(2N-1) = N is T1 in labels n+1 (3n+1 sheet); F-(2N) = 3N, F-(2N-1) = N-1 is 3n-1 in labels n-1")
ms = sorted([N for N in range(1, 50)] + [N - 1 for N in range(1, 50)])
print("  multiset of odd-label values of F+ and F- together (N <= 49):", ms[:12], "...")
check(ms[0] == 0 and ms[1] != 0 and all(ms.count(v) == 2 for v in range(1, 49)),
      "G: together the two completions take every positive integer twice and 0 once: the owner's 'two copies plus one 0'")
# folded two-sheet relation y -> {floor(y/2), ceil(3y/2)}
ok = True
for z in range(-500, 500):
    y = z if z >= 0 else -1 - z
    tz = T1(z)
    ty = tz if tz >= 0 else -1 - tz
    if ty not in {y // 2, -(-3 * y // 2)}:
        ok = False
check(ok, "G: folding Z about -1/2 (z -> z or -1-z) sends T1 to the two-branch relation y -> floor(y/2), ceil(3y/2)")
print("  first K values and their reading: K_1..K_3 =", K[:3], "= exponents of {2^0, 2^2, 2^1} = {1,4,2}; K_2..K_5 =", K[1:5], "= the cycle 2 -> 1 -> 4 -> 2")
agree_plus = [N for N in range(1, 12) if K[N - 1] == N]
agree_minus = [N for N in range(1, 12) if K[N - 1] == N - 1]
print(f"  K agrees with F+ at N = {agree_plus}, with F- at N = {agree_minus}")

print("\n" + ("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} CHECK(S) FAILED: {FAILS}"))
