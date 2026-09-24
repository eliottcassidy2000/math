#!/usr/bin/env python3
"""Owner's odd-square brackets, part 1: escape sets S_r for every multiplier, the Collatz
'microcosm', multiples k*p in p's own bracket, and tests of whether {2,3,11} is special.

Bracket B_m = ((2m-1)^2, (2m+1)^2] for m >= 1; B_0 = {1} (the degenerate bracket of 1).
S_r = {n >= 1 : n and r*n lie in the same bracket}  (r*n real; exact rational arithmetic).
rho(n) = (top of n's bracket)/n = (2m(n)+1)^2/n, so that  n in S_r  <=>  r <= rho(n).

Session collatz-procgen-20260923/24, brackets lane (procgen_brackets_20260924). Pure Python,
exact integer arithmetic; runs in about a minute.  Usage: python3 <this> [NMAX]
"""
import math, sys
from fractions import Fraction as Fr

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3 * 10**5

def m_of(n):
    """bracket index: n in B_m, i.e. (2m-1)^2 < n <= (2m+1)^2; m(1) = 0."""
    if n == 1:
        return 0
    s = math.isqrt(n - 1)                   # s^2 <= n-1 < (s+1)^2, so the smallest odd square >= n is ...
    t = s + 1 if (s + 1) % 2 == 1 else s + 2  # smallest odd t with t^2 >= n
    return (t - 1) // 2

def top(n):
    return (2 * m_of(n) + 1) ** 2

def same(a_num, a_den, n):
    """is the real number a_num/a_den in the bracket of the integer n?"""
    m = m_of(n)
    if m == 0:
        return a_num == a_den            # B_0 = {1}
    lo, hi = (2 * m - 1) ** 2, (2 * m + 1) ** 2
    return lo * a_den < a_num <= hi * a_den

def primes_upto(N):
    s = bytearray([1]) * (N + 1); s[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(N) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, N + 1, i)))
    return [i for i in range(N + 1) if s[i]]

def h_real(m):   # real-interval threshold: ((2m-1)^2, (2m+1)^2/r] nonempty  <=>  r < h(m)
    return Fr((2 * m + 1) ** 2, (2 * m - 1) ** 2)

def g_int(m):    # integer threshold: some integer n in B_m with rn in B_m  <=>  r <= g(m)
    return Fr((2 * m + 1) ** 2, (2 * m - 1) ** 2 + 1)

def mmax_formula(r):
    x = (math.sqrt(r) + 1) / (2 * (math.sqrt(r) - 1))
    return math.floor(x + 1e-12)

def M_int(r):
    m = 0
    while g_int(m + 1) >= r:
        m += 1
    return m

def S_rational(r, bound):
    p, q = r.numerator, r.denominator
    return [n for n in range(1, bound + 1) if same(p * n, q, n)]

print("=" * 100)
print("PART 1. Escape sets S_r, thresholds, the Collatz microcosm, multiples, and {2,3,11}")
print("=" * 100)

# ---------------------------------------------------------------- 1a. thresholds
print("\n(1a) Thresholds.  Real version: ((2m-1)^2,(2m+1)^2/r] nonempty iff r < h(m) = ((2m+1)/(2m-1))^2,")
print("     i.e. iff m < X(r) = (sqrt r + 1)/(2(sqrt r - 1)).  Integer version: iff r <= g(m) = (2m+1)^2/((2m-1)^2+1).")
print("     h and g are strictly decreasing in m >= 1, so S_r meets exactly the brackets m = 1..M_int(r).")
print("     m :   h(m)      g(m)")
for m in range(1, 11):
    print(f"   {m:3d}  {float(h_real(m)):8.5f}  {float(g_int(m)):8.5f}")
mono = all(h_real(m) > h_real(m + 1) and g_int(m) > g_int(m + 1) and g_int(m) < h_real(m) for m in range(1, 5000))
print("     monotone and g < h for m = 1..5000:", mono)
# compare floor formula with the exact integer threshold on a grid of rationals
over, eq_case, agree, tot = [], [], 0, 0
for qd in range(1, 41):
    for pn in range(qd + 1, 10 * qd + 1):
        r = Fr(pn, qd)
        if r.denominator != qd:
            continue
        tot += 1
        a, b = mmax_formula(float(r)), M_int(r)
        if a == b:
            agree += 1
        else:
            over.append((str(r), a, b))
        # exact-integer X(r) case r = h(m)
        for m in range(1, 60):
            if h_real(m) == r:
                eq_case.append((str(r), m))
print(f"     grid of {tot} rationals r = p/q in (1,10], q <= 40: floor-formula m_max(r) = exact M_int(r) for {agree};")
print(f"     formula overshoots by exactly 1 for {len(over)} of them (all with r in the window (g(m), h(m)] for m = m_max):")
print("       first cases (r, m_max formula, exact):", over[:12])
print("     r = h(m) exactly (X(r) an integer, real interval degenerate):", eq_case[:6])
okwin = all(g_int(b + 1) < Fr(r) <= h_real(b + 1) or Fr(r) == h_real(b + 1) for r, a, b in [(Fr(x), a, b) for x, a, b in over])
print("     every overshoot lies in (g(m_max), h(m_max)]:", okwin)
# brute-force check of M_int on the grid (S_r computed directly up to 2000 > all relevant brackets for r >= 1.1)
bad = 0
for r in [Fr(p, q) for q in range(1, 13) for p in range(q + 1, 6 * q + 1)]:
    if r <= Fr(11, 10):
        continue
    S = S_rational(r, 3000)
    mm = max([m_of(n) for n in S], default=0)
    if mm != M_int(r):
        bad += 1
print("     brute force (n <= 3000) max bracket index of S_r equals M_int(r) on all tested r > 1.1:", bad == 0)

# ---------------------------------------------------------------- 1b. tables
print("\n(1b) Escape sets S_r (all n; complete because the brute-force bound exceeds (2 M_int + 1)^2):")
for r in [Fr(2), Fr(3), Fr(3, 2), Fr(4, 3), Fr(5, 4), Fr(5, 3), Fr(4), Fr(9, 2), Fr(5)]:
    S = S_rational(r, 5000)
    by = {}
    for n in S:
        by.setdefault(m_of(n), []).append(n)
    P = [n for n in S if all(n % d for d in range(2, math.isqrt(n) + 1)) and n > 1]
    print(f"   r = {str(r):4s}  m_max formula {mmax_formula(float(r))}, exact M_int {M_int(r)}, |S_r| = {len(S)}, max {max(S) if S else '-'}")
    for m in sorted(by):
        v = by[m]
        print(f"        B_{m}: {v if len(v) <= 14 else (str(v[:6])[:-1] + ', ..., ' + str(v[-3:])[1:])}")
    print(f"        primes: {P}")

# ---------------------------------------------------------------- 1c. Collatz steps
print("\n(1c) Integer Collatz-type steps staying in their own bracket (brute force n <= NMAX; completeness")
print("     proved by the quadratic inequalities printed after each line):")
steps = {
    "plus odd  (3n+1)/2": (lambda n: n % 2 == 1, lambda n: (3 * n + 1, 2)),
    "halving   n/2     ": (lambda n: n % 2 == 0, lambda n: (n, 2)),
    "minus odd (3n-1)/2": (lambda n: n % 2 == 1, lambda n: (3 * n - 1, 2)),
    "3n+1 (unshortcut) ": (lambda n: n % 2 == 1, lambda n: (3 * n + 1, 1)),
    "3n-1 (unshortcut) ": (lambda n: n % 2 == 1, lambda n: (3 * n - 1, 1)),
    "5n+1 (unshortcut) ": (lambda n: n % 2 == 1, lambda n: (5 * n + 1, 1)),
    "(5n+1)/2 shortcut ": (lambda n: n % 2 == 1, lambda n: (5 * n + 1, 2)),
    "(5n-1)/2 shortcut ": (lambda n: n % 2 == 1, lambda n: (5 * n - 1, 2)),
    "(7n+1)/2 shortcut ": (lambda n: n % 2 == 1, lambda n: (7 * n + 1, 2)),
}
Ssteps = {}
for name, (dom, f) in steps.items():
    S = [n for n in range(1, NMAX + 1) if dom(n) and same(*f(n), n)]
    Ssteps[name] = S
    print(f"   {name}: {S}")
print("   proofs (smallest admissible n in B_m against the top (2m+1)^2):")
print("     (3n+1)/2: 3((2m-1)^2+2)+1 <= 2(2m+1)^2  <=>  m^2-5m+2 <= 0  <=>  1 <= m <= 4")
print("     (3n-1)/2: 3((2m-1)^2+2)-1 <= 2(2m+1)^2  <=>  2m^2-10m+3 <= 0 <=>  1 <= m <= 4  (plus the fixed point n = 1 in B_0)")
print("     n/2     : 2(2m-1)^2+2 <= (2m+1)^2       <=>  4m^2-12m+3 <= 0 <=>  1 <= m <= 2")
chk = all((m * m - 5 * m + 2 <= 0) == (1 <= m <= 4) and (2 * m * m - 10 * m + 3 <= 0) == (1 <= m <= 4)
          and (4 * m * m - 12 * m + 3 <= 0) == (1 <= m <= 2) for m in range(1, 10**4))
print("     inequality equivalences checked for m < 10^4:", chk)
half = Ssteps["halving   n/2     "]
print("   halving set = 2 * S_2 exactly:", half == [2 * n for n in S_rational(Fr(2), 5000)])
T_all = sorted(set(Ssteps["plus odd  (3n+1)/2"]) | set(half))
U_all = sorted(set(Ssteps["minus odd (3n-1)/2"]) | set(half))
print(f"   MICROCOSM.  plus sheet T: n with T(n) in n's bracket = {T_all}; largest = {max(T_all)}.")
print(f"               minus sheet : {U_all}; largest = {max(U_all)}.")
print("               => every move of either shortcut map from n >= 54 (in particular every n > 81 = 9^2,")
print("                  brackets m >= 5) changes bracket (PROVED above).")
# bracket-index change per step for large n
import random
random.seed(1)
rat_up, rat_dn = [], []
for _ in range(20000):
    n = random.randrange(10**8, 10**9)
    if n % 2:
        rat_up.append(m_of((3 * n + 1) // 2) / m_of(n))
    else:
        rat_dn.append(m_of(n // 2) / m_of(n))
print(f"   large n (1e8..1e9): bracket index ratio after an odd step {sum(rat_up)/len(rat_up):.5f} (sqrt(3/2) = {math.sqrt(1.5):.5f}),"
      f" after halving {sum(rat_dn)/len(rat_dn):.5f} (1/sqrt 2 = {1/math.sqrt(2):.5f})")
# pairs {2i-1,2i} both of whose moves stay inside B_1..: the fully internal pairs
full_pairs = [i for i in range(1, 100) if (2 * i - 1) in T_all and (2 * i) in T_all]
print(f"   pairs {{2i-1,2i}} whose up- AND down-move both stay in the bracket (plus sheet): i in {full_pairs}"
      f"  (pairs {{3,4}}, {{5,6}} inside B_1)")
full_pairs_minus = [i for i in range(1, 100) if (2 * i + 1) in U_all and (2 * i) in U_all]
print(f"   same for the minus-sheet pairs {{2i,2i+1}}: i in {full_pairs_minus}")

# ---------------------------------------------------------------- 1d. rho and multiples
print("\n(1d) The single function rho(n) = (2m(n)+1)^2/n encodes every S_r:  S_r = {n : rho(n) >= r}.")
rho = sorted(((Fr(top(n), n), n) for n in range(2, 400)), reverse=True)
print("     escape order (n by decreasing rho), first 24:", [(n, round(float(x), 4)) for x, n in rho[:24]])
PR = primes_upto(10**6)
prho = sorted(((Fr(top(p), p), p) for p in PR if p < 2000), reverse=True)
print("     primes by decreasing rho, first 16:", [(p, round(float(x), 4)) for x, p in prho[:16]])
print("     => prime escape sets P_r: r in (3, 4.5] -> {2};  (2.2727, 3] -> {2,3};  (1.9231, 2.2727] -> {2,3,11};")
print("        then 13 (1.923), 5 (1.8), 29 (1.690), 31 (1.581), 53 (1.528), 17 (1.471), 83 (1.458), ...")
kp = [(k, p) for p in PR[:2000] for k in range(2, 10) if same(k * p, 1, p)]
print("     all (k, p), k >= 2 integer, p prime < 10^6 (search k <= 9 suffices: k <= rho(p) <= 4.5), with kp in p's bracket:", kp)
kn = [(k, n) for n in range(2, 10**5) for k in range(2, 6) if same(k * n, 1, n)]
print("     same for all integers n < 10^5:", kn)
print("     singleton reading: every prime p outside {2,3,11} is the ONLY multiple of p in its own bracket;")
print("     for 2, 3, 11 the bracket also holds 2p (and 3p for p = 2, 3; 4p for p = 2).")

# multiples of p across later brackets: skips
def mult_count(p, j):
    return (2 * j + 1) ** 2 // p - (2 * j - 1) ** 2 // p

skipfree, first_skip = [], {}
for p in PR:
    if p > 20000:
        break
    m0 = m_of(p)
    j = m0
    skipped = None
    while 8 * j < p + 8:                      # once 8j >= p every bracket has a multiple
        if j > m0 and mult_count(p, j) == 0:
            skipped = j
            break
        j += 1
    if skipped is None:
        skipfree.append(p)
    else:
        first_skip[p] = skipped
print("     skip-free primes (every bracket B_j, j >= m(p), contains a multiple of p), p < 20000:", skipfree)
print("       count", len(skipfree), "; largest", max(skipfree), "(e.g. 41 skips: 41 in B_3, 82 in B_5, B_4 = (49,81] holds no multiple of 41)")
print("       PROOF of finiteness: if p in B_m then all multiples of p up to (4m+1)^2 lie in B_m..B_2m and number at most")
print("       (4m+1)^2/p < (4m+1)^2/(2m-1)^2, which is < m+1 (the number of those brackets) for every m >= 5; so every")
print("       prime above 81 skips a bracket, and the list above (from the finite check p <= 81) is complete.")
print("       m >= 5 criterion checked:", all((4*m+1)**2 < (m+1)*(2*m-1)**2 for m in range(5, 10**4)), "; fails at m = 4:", (4*4+1)**2 < 5*7**2)
# the sieve-window view
print("\n(1e) Sieve-window view: for m >= 2 the primes of B_m are exactly the n in B_m coprime to all primes <= 2m+1,")
print("     and the new sieving prime q = 2m+1 (if prime) first removes q^2 = the top of B_m.  Mertens: the naive")
print("     sieve prediction 8m*prod_{p<=2m+1}(1-1/p) overshoots the true count by 2e^-gamma = 1.1229 (classical).")
Mb = 4000
LIM = (2 * Mb + 1) ** 2
sieve = bytearray([1]) * (LIM + 1); sieve[0:2] = b"\x00\x00"
for i in range(2, math.isqrt(LIM) + 1):
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(range(i * i, LIM + 1, i)))
prod, ps, idx = 1.0, [p for p in PR if p <= 2 * Mb + 1], 0
ratios = []
for m in range(1, Mb + 1):
    while idx < len(ps) and ps[idx] <= 2 * m + 1:
        prod *= 1 - 1 / ps[idx]; idx += 1
    lo, hi = (2 * m - 1) ** 2 + 1, (2 * m + 1) ** 2
    c = sum(sieve[lo:hi + 1])
    if m >= 100:
        ratios.append(c / (8 * m * prod))
gam = 0.5772156649015329
print(f"     mean over 100 <= m <= {Mb} of  #primes(B_m) / (8m prod(1-1/p)) = {sum(ratios)/len(ratios):.4f};"
      f"  e^gamma/2 = {math.exp(gam)/2:.4f}")
for lo_, hi_ in [(100, 500), (500, 1500), (1500, 4000)]:
    rr = ratios[lo_ - 100:hi_ - 100]
    print(f"       m in [{lo_},{hi_}): {sum(rr)/len(rr):.4f}")

# ---------------------------------------------------------------- 1f. is {2,3,11} special?
print("\n(1f) Is anything special about {2,3,11} beyond m <= 2?")
wief = [p for p in primes_upto(2 * 10**6) if p != 3 and pow(3, p - 1, p * p) == 1]
print("     base-3 Wieferich (Mirimanoff) primes p < 2*10^6:", wief, " (3^5 = 243 = 2*11^2 + 1)")
sq = [k for k in range(1, 400) if math.isqrt((3**k - 1) // 2) ** 2 == (3**k - 1) // 2]
print("     k <= 400 with (3^k - 1)/2 a square:", sq, "-> values", [(3**k - 1) // 2 for k in sq],
      " (Ljunggren 1943: these are all, for every k)")
print("     E-graph climb from 1 (x -> 3x+1): 1, 4, 13, 40, 121, 364, ... = (3^k-1)/2; its odd squares are 1 and 121 = 11^2 only.")
brk = [m for m in range(1, 10**5) if math.isqrt((2 * m - 1) ** 2) ** 2 == (2 * m - 1) ** 2
       and any((2 * m - 1) ** 2 == 3**a and (2 * m + 1) ** 2 == (3**(a + 1) - 1) // 2 for a in range(1, 25))]
print("     brackets B_m = (3^a, (3^(a+1)-1)/2] for m < 10^5:", brk, "(B_5 = (81,121]; unique for all m by Ljunggren)")
trunk_minus = [(2**(2 * n + 1) + 1) // 3 for n in range(0, 12)]
trunk_plus = [(4**n - 1) // 3 for n in range(1, 12)]
isp = lambda x: x > 1 and all(x % d for d in range(2, math.isqrt(x) + 1))
print("     minus trunk (2^(2n+1)+1)/3:", trunk_minus[:8], "; its primes (Wagstaff):", [t for t in trunk_minus if isp(t)])
print("     trunk members in S_2:", [t for t in trunk_minus + trunk_plus if t in S_rational(Fr(2), 100)],
      "; all odd members of S_2:", [n for n in S_rational(Fr(2), 100) if n % 2])
print("     => S_2's odd part {3,11} equals the minus trunk inside (1, 12.5] because (1,12.5] holds only the")
print("        trunk values 1, 3, 11; a random 2-subset of the odd primes {3,5,7,11} hits {3,11} with chance 1/6.")
# every small prime is Wieferich to some small base: smallest base b >= 2, p not dividing b, with b^(p-1) = 1 mod p^2
sm = []
for p in primes_upto(60):
    b = next(b for b in range(2, p * p + 2) if b % p != 0 and pow(b, p - 1, p * p) == 1)
    sm.append((p, b))
print("     smallest base b >= 2 (p not dividing b) with b^(p-1) = 1 mod p^2, p < 60:", sm)
print("     the solutions mod p^2 are the p-1 roots of unity, one in each nonzero class mod p, so a fixed base is")
print("     Wieferich for p with chance ~ 1/p: '11 is Wieferich to base 3' is a 1-in-11 event, and 2 (base 5) and")
print("     3 (base 8) are Wieferich to small bases too.  No mechanism ties either fact to 2p <= (2m+1)^2.")
print("\n(1g) Figurate identities at the bracket tops (exact, elementary):")
ok1 = all(((3 * (2 * m + 1)**2 + 1) // 2) // 2 == (m + 1)**3 - m**3 for m in range(1, 10**4))
ok2 = all((3 * (2 * m + 1)**2 - 1) // 2 == 6 * (m + 1) * m + 1 for m in range(1, 10**4))
print("     T^2((2m+1)^2) = (3(2m+1)^2+1)/4 = 3m^2+3m+1 = (m+1)^3 - m^3 (centred hexagonal), m < 10^4:", ok1)
print("     U((2m+1)^2) = (3(2m+1)^2-1)/2 = 6m(m+1)+1 (star number), m < 10^4:", ok2)
print("     so every odd square has stopping time 2 under T and lands in bracket ~ (sqrt3/2) m.")
