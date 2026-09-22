#!/usr/bin/env python3
"""Lane scaffolding_audit (session collatz-mod6-20260917, mac-mini; written 2026-09-21).

Exact audit of every claim in the pasted "peripheral discoveries" block (the
scaffolding of the geometric Collatz blueprint refuted as a proof in
05-knowledge/results/collatz_blueprint_20260921_synthesis.md).

Probes P1..P10 follow the numbered claims of the paste.  Every check uses
explicit `raise` so it survives `python -O`.  Exact integer / Fraction
arithmetic throughout; sympy is used only for factorisation (P6, P7) and the
results are re-checked by multiplication.  RAM << 1 GB, runtime ~1 minute.

Reproduce:
  python3 04-computation/experiments/collatz_mod6_20260917_scaffolding_audit.py \
      > 05-knowledge/results/collatz_mod6_20260917_scaffolding_audit.out
"""
import itertools
import math
from fractions import Fraction

import sympy

FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        raise AssertionError(msg)


def T(n):
    return n * (n + 1) // 2


def binom(n, k):
    # C(n,0)=1 for every integer n (so S_0(n)=C(n-1,0)=1, including n=0); otherwise ordinary.
    if k == 0:
        return 1
    if k < 0 or n < 0 or k > n:
        return 0
    return math.comb(n, k)


def S(d, n):
    """d-simplex number S_d(n) = C(n+d-1, d): S_1=n, S_2=T, S_3=Te, S_4=Pt."""
    return binom(n + d - 1, d)


def Te(n):
    return S(3, n)


def Pt(n):
    return S(4, n)


# ---------------------------------------------------------------- polynomials
# exact bivariate polynomials in (A,B): dict {(i,j): Fraction}
def padd(p, q, sign=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + sign * v
    return {k: v for k, v in r.items() if v != 0}


def pmul(p, q):
    r = {}
    for (i1, j1), v1 in p.items():
        for (i2, j2), v2 in q.items():
            k = (i1 + i2, j1 + j2)
            r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}


def pconst(c):
    return {(0, 0): Fraction(c)} if c != 0 else {}


PA = {(1, 0): Fraction(1)}
PB = {(0, 1): Fraction(1)}


def ptri(p):
    """T(p) = p(p+1)/2 as a polynomial."""
    return pmul(pmul(p, padd(p, pconst(1))), pconst(Fraction(1, 2)))


def pshift(p, c):
    return padd(p, pconst(c))


def pstr(p):
    if not p:
        return "0"
    terms = []
    for (i, j) in sorted(p, key=lambda k: (-(k[0] + k[1]), -k[0])):
        v = p[(i, j)]
        mono = ("A^%d" % i if i > 1 else ("A" if i == 1 else "")) + ("B^%d" % j if j > 1 else ("B" if j == 1 else ""))
        if mono == "":
            terms.append(str(v))
        elif v == 1:
            terms.append(mono)
        elif v == -1:
            terms.append("-" + mono)
        else:
            terms.append("%s*%s" % (v, mono))
    return " + ".join(terms).replace("+ -", "- ")


def hdr(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


print("scaffolding_audit: exact audit of the pasted 'peripheral discoveries' block")
print("session collatz-mod6-20260917 (mac-mini), lane written 2026-09-21")

# =============================================================================
hdr("P1. Claim 1: 'core local identity' T((A+1)(B+1)-1)-T(AB-1) = T(A)+T(B)+Shear(A,B),"
    " Shear = A(B^2-1)+B(A^2-1); horizontal/vertical flows of S(A,B)")


def shear(A, B):
    return A * (B * B - 1) + B * (A * A - 1)


fails = 0
for A in range(1, 8):
    for B in range(1, 8):
        lhs = T((A + 1) * (B + 1) - 1) - T(A * B - 1)
        rhs = T(A) + T(B) + shear(A, B)
        if lhs != rhs:
            fails += 1
print("(i) stated identity, (A,B) in [1,7]^2: %d of 49 pairs FAIL" % fails)
check(fails == 49, "expected the stated identity to fail on all 49 pairs")
print("    witness A=B=1: LHS = T(3)-T(0) = %d, RHS = %d" % (T(3) - T(0), T(1) + T(1) + shear(1, 1)))
check(T(3) - T(0) == 6 and T(1) + T(1) + shear(1, 1) == 2, "witness A=B=1 must be 6 vs 2")

# exact polynomial identity for the left side
N = padd(pmul(PA, PB), padd(PA, PB))               # N = AB+A+B = (A+1)(B+1)-1
LHS = padd(ptri(N), ptri(pshift(pmul(PA, PB), -1)), -1)
RHS_paste = padd(padd(ptri(PA), ptri(PB)), padd(pmul(PA, padd(pmul(PB, PB), pconst(-1))),
                                                 pmul(PB, padd(pmul(PA, PA), pconst(-1)))))
diff = padd(LHS, RHS_paste, -1)
print("(ii) exact polynomials:")
print("    LHS  = T(AB+A+B)-T(AB-1) = %s" % pstr(LHS))
print("    RHS  = T(A)+T(B)+A(B^2-1)+B(A^2-1) = %s" % pstr(RHS_paste))
print("    LHS-RHS = %s  = A(B+1)+B(A+1)" % pstr(diff))
check(diff == padd(pmul(PA, pshift(PB, 1)), pmul(PB, pshift(PA, 1))), "LHS-RHS must equal A(B+1)+B(A+1)")
# closed forms for the left side
half = pconst(Fraction(1, 2))
closed = pmul(pmul(pshift(padd(PA, PB), 1), padd(pmul(pconst(2), pmul(PA, PB)), padd(PA, PB))), half)
check(LHS == closed, "LHS must equal (A+B+1)(2AB+A+B)/2")
true_shear = pmul(pmul(PA, PB), pshift(padd(PA, PB), 2))          # AB(A+B+2)
check(LHS == padd(padd(ptri(PA), ptri(PB)), true_shear), "LHS must equal T(A)+T(B)+AB(A+B+2)")
check(true_shear == padd(pmul(pconst(2), pmul(PA, ptri(PB))), pmul(pconst(2), pmul(PB, ptri(PA)))),
      "AB(A+B+2) must equal 2A*T(B)+2B*T(A)")
print("    PROVED: T((A+1)(B+1)-1)-T(AB-1) = sum_{k=AB}^{AB+A+B} k = (A+B+1)(2AB+A+B)/2")
print("                                  = T(A)+T(B)+AB(A+B+2) = T(A)+T(B)+2A*T(B)+2B*T(A).")

# (iii) which shifted readings make the stated FORM true?  LHS' = T(N+s1)-T(AB+s2), RHS' = T(A+t1)+T(B+t2)+Shear
sols = []
for s1, s2, t1, t2 in itertools.product(range(-3, 4), repeat=4):
    ok = True
    for A in range(0, 8):
        for B in range(0, 8):
            if T((A + 1) * (B + 1) - 1 + s1) - T(A * B + s2) != T(A + t1) + T(B + t2) + shear(A, B):
                ok = False
                break
        if not ok:
            break
    if ok:
        sols.append((s1, s2, t1, t2))
print("(iii) shifted readings T(N+s1)-T(AB+s2) = T(A+t1)+T(B+t2)+Shear valid on [0,7]^2, shifts in [-3,3]:")
print("    solutions (s1,s2,t1,t2): %s" % sols)
check(sols == [(-3, -2, -2, -2), (-1, 0, 0, 0)], "exactly two shifted readings expected")
# polynomial proof of the repaired identity
LHS_rep = padd(ptri(pshift(N, -1)), ptri(pmul(PA, PB)), -1)
check(LHS_rep == RHS_paste, "repaired identity T((A+1)(B+1)-2)-T(AB) = T(A)+T(B)+Shear must hold as polynomials")
print("    PROVED (polynomial identity): T((A+1)(B+1)-2) - T(AB) = T(A)+T(B)+A(B^2-1)+B(A^2-1)")
print("      i.e. the paste is off by one in BOTH triangular arguments; the second solution is its")
print("      reflection under T(n-1)=T(-n).  Witness of the repair: A=B=1: T(2)-T(1)=%d = T(1)+T(1)+0." % (T(2) - T(1)))
# general lemma: no quadratic reading T(n)=a n^2+b n+c makes the stated form true (AB-coefficient 2a vs 0 with a=1/2 forced)
print("    PROVED (coefficient comparison): for ANY quadratic T(n)=a n^2+b n+c the stated form forces a=1/2")
print("      from the A^2B term and then leaves an uncancelled 2AB on the left; so no quadratic reading works.")

# (iv) flows
def g_h(A, B):
    return A + B * B + 2 * A * B + B


def g_v(A, B):
    return B + A * A + 2 * A * B + A


# compatibility
comp_ok = all(g_h(A, B + 1) - g_h(A, B) == g_v(A + 1, B) - g_v(A, B) for A in range(0, 12) for B in range(0, 12))
check(comp_ok, "flows must be compatible")
S_poly = padd(RHS_paste, {})   # candidate S = T(A)+T(B)+Shear = AB(A+B)+C(A,2)+C(B,2)


def S_val(A, B):
    return A * B * (A + B) + binom(A, 2) + binom(B, 2)


flow_ok = all(S_val(A + 1, B) - S_val(A, B) == g_h(A, B) and S_val(A, B + 1) - S_val(A, B) == g_v(A, B)
              for A in range(0, 12) for B in range(0, 12))
check(flow_ok, "S = AB(A+B)+C(A,2)+C(B,2) must satisfy both flows")
check(S_poly == padd(pmul(pmul(PA, PB), padd(PA, PB)), padd(ptri(pshift(PA, -1)), ptri(pshift(PB, -1)))),
      "T(A)+T(B)+Shear must equal AB(A+B)+C(A,2)+C(B,2)")
# LHS's own flows
lhs_h = [(T((A + 2) * (B + 1) - 1) - T((A + 1) * B - 1)) - (T((A + 1) * (B + 1) - 1) - T(A * B - 1)) - g_h(A, B)
         for A in range(0, 6) for B in range(0, 6)]
lhs_h_extra = set((A, B, (T((A + 2) * (B + 1) - 1) - T((A + 1) * B - 1)) - (T((A + 1) * (B + 1) - 1) - T(A * B - 1)) - g_h(A, B) - (2 * B + 1))
                  for A in range(0, 6) for B in range(0, 6))
check(all(x[2] == 0 for x in lhs_h_extra), "LHS horizontal flow must be g_h + 2B+1")
print("(iv) flows Delta_A S = A+B^2+2AB+B, Delta_B S = B+A^2+2AB+A are compatible (mixed differences both 2A+2B+2)")
print("    and their UNIQUE solution up to the constant S(0,0) is")
print("      S(A,B) = AB(A+B) + C(A,2) + C(B,2) = T(A)+T(B)+A(B^2-1)+B(A^2-1)  (the paste's RIGHT side),")
print("      = T((A+1)(B+1)-2) - T(AB).  The paste's LEFT side T((A+1)(B+1)-1)-T(AB-1) has horizontal flow")
print("      Delta_A = (A+B^2+2AB+B) + (2B+1), so it is NOT the S of the flows.  Verdict: TRUE-with-repair.")
print("    uniqueness: two solutions of both difference equations differ by a function constant in A and in B.")

# =============================================================================
hdr("P2. Claim 2: tetrahedral and pentatope laws; general d-simplex convolution law")
bad3 = [(A, B) for A in range(0, 12) for B in range(0, 12)
        if Te(A + B + 1) != Te(A) + Te(B) + (A + 1) * T(B) + (B + 1) * T(A) + (A + 1) * (B + 1)]
bad4 = [(A, B) for A in range(0, 12) for B in range(0, 12)
        if Pt(A + B + 1) != Pt(A) + Pt(B) + (A + 1) * Te(B) + (B + 1) * Te(A) + T(A + 1) * T(B + 1)]
check(bad3 == [] and bad4 == [], "tetrahedral / pentatope laws must hold on [0,11]^2")
print("(i) FINITE-EXACT [0,11]^2: tetrahedral law holds (%d failures), pentatope law holds (%d failures)." % (len(bad3), len(bad4)))
# 2D law for comparison
bad2 = [(A, B) for A in range(0, 12) for B in range(0, 12) if T(A + B + 1) != T(A) + T(B) + (A + 1) * (B + 1)]
check(bad2 == [], "2D law T(A+B+1)=T(A)+T(B)+(A+1)(B+1)")
print("    the 2D member of the same family, T(A+B+1) = T(A)+T(B)+(A+1)(B+1), also holds (this is the")
print("    correct 'locking' identity that Claim 1 garbled).")
# Vandermonde proof, checked exactly for all d<=8, all splits a+b=d, on [0,11]^2
for d in range(1, 9):
    for a in range(0, d + 1):
        b = d - a
        for A in range(0, 12):
            for B in range(0, 12):
                lhs = binom(A + B + d, d)
                rhs = sum(binom(A + a, i) * binom(B + b, d - i) for i in range(0, d + 1))
                check(lhs == rhs, "Vandermonde split failed d=%d a=%d" % (d, a))
print("(ii) PROVED (Chu-Vandermonde): for every d>=1 and every split a+b=d,")
print("      S_d(A+B+1) = C(A+B+d,d) = sum_{i=0}^{d} C(A+a,i) C(B+b,d-i),  S_d(n):=C(n+d-1,d).")
print("    Proof: choose a d-subset of a set of A+B+d = (A+a)+(B+b) elements by how many lie in the first block.")
print("    Tetrahedral law = split (a,b)=(2,1) rewritten with Pascal: C(B+1,3)+(A+2)C(B+1,2) = Te(B)+(A+1)T(B),")
print("      C(A+2,2)(B+1) = (B+1)T(A)+(A+1)(B+1), C(A+2,3) = Te(A).")
print("    Pentatope law = split (2,2): C(B+2,4)+(A+2)C(B+2,3) = Pt(B)+(A+1)Te(B) (and symmetrically),")
print("      middle term C(A+2,2)C(B+2,2) = T(A+1)T(B+1).")
# symmetric general simplex convolution: S_d(A+B+1) = sum_{i+j=d} S_i(A+1) S_j(B)
for d in range(1, 9):
    for A in range(0, 12):
        for B in range(0, 12):
            check(S(d, A + B + 1) == sum(S(i, A + 1) * S(d - i, B) for i in range(0, d + 1)),
                  "simplex convolution failed d=%d" % d)
print("(iii) PROVED (the identity sum_k C(x+k,k)C(y+d-k,d-k) = C(x+y+d+1,d), x=A, y=B-1): the simplex convolution")
print("      S_d(A+B+1) = sum_{i+j=d} S_i(A+1) S_j(B), checked exactly for d<=8 on [0,11]^2.")
print("    Physical content that survives: a convolution (set-partition count) identity; there is no shear,")
print("    rotation, or 'locking plate' -- every term is C(x,i)C(y,d-i) for a two-block split.")

# =============================================================================
hdr("P3. Claim 3: the g-operator")
print("(i) 'AgB = AB*(AgB)': as an equation in the unknown AgB it forces (AB-1)*(AgB)=0, so AgB=0 unless AB=1:")
print("    circular / vacuous.  REFUTED as a definition.")
incons = [(B, B ** B, math.factorial(B)) for B in range(1, 7) if B ** B != math.factorial(B)]
print("(ii) '(AB)^(AB)' at A=1 gives B^B; '1gB=B!' gives B!; B^B != B! for B>=2: witnesses %s" % incons[:3])
check(incons[0] == (2, 4, 2), "B=2 witness 4 vs 2")
print("    REFUTED (mutually inconsistent readings; minimal witness B=2: 4 vs 2).")
# consistent readings
def g_fact(A, B):
    return math.factorial(A) * math.factorial(B)


ok_rec = all(g_fact(A, B) == A * B * g_fact(A - 1, B - 1) for A in range(1, 8) for B in range(1, 8))
ok_1 = all(g_fact(1, B) == math.factorial(B) for B in range(0, 8))
ok_sym = all(g_fact(A, A) == math.factorial(A) ** 2 for A in range(0, 8))
check(ok_rec and ok_1 and ok_sym, "A!B! reading must satisfy the recursion and both boundary claims")
print("(iii) PROVED: the reading AgB := A!*B! satisfies the (presumably intended) recursion AgB = AB*((A-1)g(B-1)),")
print("      1gB = B!, and AgA = (A!)^2 simultaneously.  Any solution f of f(A,B)=AB f(A-1,B-1) with f(1,B)=B!")
print("      has f(0,k)=k! and f(A,A)=(A!)^2 f(0,0)=(A!)^2, so A!B! is the unique symmetric solution.")
print("      The alternative AgB := A^B B! (recursion AgB = AB*(Ag(B-1))) gives 1gB=B! but AgA = A^A A! != (A!)^2")
print("      (A=2: %d vs %d).  '(AB)^(AB)' satisfies neither recursion." % (2 ** 2 * 2, 4))
print("(iv) 'squares the information capacity, mapping the Pell and square-triangular families': SCOPE, no map found")
print("      ((A!)^2 is a square by construction; square-triangular numbers T(m)=k^2 are the Pell equation")
print("      (2m+1)^2-8k^2=1 and involve no factorials).")

# =============================================================================
hdr("P4. Claim 4: tournament 'compression' T_{n-2}; Redei; counting")
for n in range(2, 13):
    check(binom(n, 2) - (n - 1) == T(n - 2) == binom(n - 1, 2), "C(n,2)-(n-1)=T(n-2)")
print("(i) TRUE (arithmetic): C(n,2)-(n-1) = C(n-1,2) = T(n-2) for all n>=2 (checked n<=12; identity C(n,2)=C(n-1,2)+(n-1)).")


def ham_paths(n, adj):
    """number of directed Hamiltonian paths of the tournament with adjacency bitmask adj[v] (out-neighbours)."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c == 0:
                continue
            nxt = adj[v] & ~mask
            while nxt:
                w = (nxt & -nxt).bit_length() - 1
                nxt &= nxt - 1
                dp[mask | (1 << w)][w] += c
    return sum(dp[full])


def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    for bits in range(1 << m):
        adj = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
        yield adj


print("(ii) labelled tournaments, exact census:")
print("      n  #T=2^C(n,2)  sum_T h(T)  n!*2^T(n-2)  min h  max h  all h odd  h-values seen")
for n in range(2, 7):
    hs = []
    for adj in all_tournaments(n):
        hs.append(ham_paths(n, adj))
    tot = sum(hs)
    expect = math.factorial(n) * 2 ** T(n - 2)
    check(tot == expect, "sum of h must be n! 2^T(n-2) at n=%d" % n)
    check(all(h % 2 == 1 for h in hs), "Redei: h odd, n=%d" % n)
    check(len(hs) == 2 ** binom(n, 2), "count of tournaments")
    seen = sorted(set(hs))
    check(7 not in seen and 21 not in seen, "h-spectrum holes 7,21 (THM-1745) n=%d" % n)
    print("     %2d  %10d  %10d  %11d  %5d  %5d  %s  %s" % (n, len(hs), tot, expect, min(hs), max(hs), all(h % 2 for h in hs), seen))
print("    PROVED (double counting): every ordering pi of [n] is the Hamiltonian path of exactly 2^{T(n-2)} labelled")
print("      tournaments (the n-1 path edges are forced, the other T(n-2) are free), so sum_T h(T) = n! 2^{T(n-2)}")
print("      and the average number of Hamiltonian paths is n!/2^{n-1} (n=6: %s)." % Fraction(math.factorial(6), 2 ** 5))
print("(iii) REFUTED ('path edges carry 0 bits'): the code (pi, T(n-2) free bits) has n!*2^{T(n-2)} codewords for")
print("      2^{C(n,2)} tournaments; it is surjective and (n!/2^{n-1})-to-one on average, so it costs")
print("      log2(n!) + T(n-2) = C(n,2) + log2(n!/2^{n-1}) >= C(n,2) bits, with strict excess for n>=3.")
for n in (3, 4, 5, 6, 8, 10):
    excess = math.log2(math.factorial(n)) - (n - 1)
    print("        n=%2d: C(n,2)=%3d, T(n-2)=%3d, log2(n!)=%.3f, excess over C(n,2) = %.3f bits" % (n, binom(n, 2), T(n - 2), math.log2(math.factorial(n)), excess))
    check(excess > 0, "excess positive n>=3")
print("      The uniform distribution on labelled tournaments has entropy exactly C(n,2) bits; no lossless code")
print("      beats it (Shannon source coding, CITED: Shannon 1948).  What IS true: h is odd (Redei 1934, CITED),")
print("      h-spectrum = odds minus {7,21} and |Aut(T)| divides h (01-canon/theorems/THM-1745-leaf-graded-...).")
print("(iv) 'defects at diffraction peaks of the quasicrystalline Wythoff wave', 'topological genus of randomness':")
print("      SCOPE -- no defined object, no map found.")

# =============================================================================
hdr("P5. Claim 5a: horizons and the '63 bottleneck' (Bang / Catalan)")
print("(i) 'horizon lines' (Z,2Z),(Z,Z^2),(W,W^3),(W,W^4) are the diagonals X=Y of X+Y, X*Y, X*Y*W... : definitions only (SCOPE).")
check(63 == 3 * 7 * 3 == 3 ** 2 * 7 == 2 ** 6 - 1 == (2 ** 3 - 1) * (2 ** 3 + 1), "63 arithmetic")
print("(ii) 63 = 3*7*3 = 3^2*7 = 2^6-1 = (2^3-1)(2^3+1) = 7*9; 9 = 3^2 = 2^3+1 (Catalan 3^2-2^3=1, Mihailescu 2004 CITED).")
# primitive prime divisors of 2^n-1, n<=40 (re-derivation of the zsigmondy lane's table, same universe)
def primitive_primes(a, n):
    val = a ** n - 1
    fac = sympy.factorint(val)
    prod = 1
    for p, e in fac.items():
        prod *= p ** e
    check(prod == val, "factorisation check")
    prim = [p for p in fac if sympy.n_order(a, p) == n]
    return fac, prim


exc2 = []
for n in range(1, 41):
    fac, prim = primitive_primes(2, n)
    if not prim:
        exc2.append(n)
print("    2^n-1, 1<=n<=40: indices with NO primitive prime divisor = %s" % exc2)
check(exc2 == [1, 6], "Bang exceptions must be n=1 (2^1-1=1) and n=6")
print("    PROVED (Bang 1886 / Zsigmondy 1892, CITED): the ONLY n>1 for which 2^n-1 lacks a primitive prime is n=6, i.e. 63.")
print("    So the paste's '63 short-circuits Bang's theorem' is TRUE-with-repair: the mechanism is 2^3+1 = 3^2 (a prime")
print("    power already dividing 2^2-1), not a 'loop switching a multiplication edge and an addition edge'; no such")
print("    loop is defined.  Inherited: wave-one zsigmondy lane, 05-knowledge/results/collatz_mod6_20260917_zsigmondy_triad.out")
print("    (Catalan reading of n=6; no prime p has ord_p(2)=6).")
check(sympy.n_order(2, 9) == 6, "ord_9(2)=6")
print("    ord_9(2) = 6 and 9 | 63; no prime has multiplicative order 6 to base 2 (primes of order 6 would be primitive).")

# =============================================================================
hdr("P6. Claim 5b: the '341 bottleneck' (Zsigmondy for 4^n-1; Cipolla pseudoprimes; the R-trunk)")
check(341 == 11 * 31 == (4 ** 5 - 1) // 3 == sum(4 ** i for i in range(5)), "341 arithmetic")
print("(i) 341 = 11*31 = (4^5-1)/3 = 4^4+4^3+4^2+4+1 (base-4 repunit): TRUE arithmetic.")
fac, prim = primitive_primes(4, 5)
print("    4^5-1 = 1023 = %s; primitive primes (ord_p(4)=5): %s" % (dict(fac), prim))
check(sorted(prim) == [11, 31], "341 has TWO primitive primes")
exc4 = []
print("    4^n-1, 1<=n<=40: primitive prime divisors")
for n in range(1, 41):
    fac, prim = primitive_primes(4, n)
    if not prim:
        exc4.append(n)
    if n <= 12 or n in (20, 30, 40):
        print("      n=%2d  4^n-1 = %s  primitive: %s" % (n, sympy.factorint(4 ** n - 1), prim))
print("    indices without a primitive prime: %s" % exc4)
check(exc4 == [], "4^n-1 has a primitive prime for every n<=40")
print("    PROVED (Zsigmondy, CITED; exceptions are only n=1 with a-b=1, n=2 with a+b a power of 2, and (2,1,6)):")
print("      4-1=3 != 1 and 4+1=5 is not a power of 2, so 4^n-1 has a primitive prime divisor for EVERY n>=1.")
print("    REFUTED: '341 stalls primitive prime generation of the 4x+1 sequence' -- minimal witness: 4^5-1 has")
print("      two primitive primes 11 and 31, both dividing 341.")
# Cipolla
print("(ii) Cipolla (1904, CITED): for prime p>=5, N_p=(4^p-1)/3 is a composite base-2 Fermat pseudoprime.")
print("      p   N_p            composite  2^(N-1) mod N  factorisation")
for p in [5, 7, 11, 13, 17, 19, 23]:
    Np = (4 ** p - 1) // 3
    comp = not sympy.isprime(Np)
    fer = pow(2, Np - 1, Np)
    check(comp and fer == 1, "Cipolla at p=%d" % p)
    print("     %3d  %-14d  %s  %d  %s" % (p, Np, comp, fer, sympy.factorint(Np)))
check((4 ** 5 - 1) // 3 == 341 and (4 ** 7 - 1) // 3 == 5461 and (4 ** 11 - 1) // 3 == 1398101, "341, 5461, 1398101")
print("    Proof sketch: N=(2^p-1)(2^p+1)/3 is composite for p>=5 (both factors >3, and 3 divides exactly one of them);")
print("      N-1 = 4(4^{p-1}-1)/3 is divisible by 2p (p | 4^{p-1}-1 by Fermat, and 2 | it), and ord_N(2) divides 2p")
print("      since 2^{2p} = 4^p = 3N+1 = 1 mod N; hence 2^{N-1} = 1 mod N.")
# R-trunk
R = lambda n: 4 * n + 1
orb = [0]
for _ in range(6):
    orb.append(R(orb[-1]))
print("(iii) R(n)=4n+1 orbit of 0: %s = (4^j-1)/3, j=0..6" % orb)
check(orb == [(4 ** j - 1) // 3 for j in range(7)], "R-orbit is (4^j-1)/3")
check(all(3 * x + 1 == 4 ** j for j, x in enumerate(orb) if j >= 1), "3n+1 of trunk = 4^j")
print("      3*((4^j-1)/3)+1 = 4^j: these are exactly the odd predecessors of the powers of two (the trunk of 1),")
print("      the inverse-fibre braid n -> R(n) of 05-knowledge/results/arithmetic_braids_20260917_collatz.md (B1)-(B2).")
print("      341 = R^5(0) = R^4(1) is on the trunk; 'zero-friction state in G_4' has no defined content (SCOPE).")

# =============================================================================
hdr("P7. Claim 5c: rational 3-cycle x^2-29/16 at x0=-7/4")
f = lambda x: x * x - Fraction(29, 16)
x0 = Fraction(-7, 4)
cyc = [x0]
for _ in range(3):
    cyc.append(f(cyc[-1]))
print("(i) orbit: %s" % [str(c) for c in cyc])
check(cyc[3] == x0 and cyc[1] == Fraction(5, 4) and cyc[2] == Fraction(-1, 4), "3-cycle -7/4 -> 5/4 -> -1/4 -> -7/4")
pts = sorted(cyc[:3])
check(pts[1] - pts[0] == pts[2] - pts[1] == Fraction(3, 2), "cycle points form an AP with difference 3/2")
print("    FINITE-EXACT: {-7/4, -1/4, 5/4} is an exact 3-cycle and an arithmetic progression with difference 3/2.")
print("    What is PROVED in the repo (cite, do not re-prove): 01-canon/theorems/THM-4139-*: the complete rational")
print("      preperiodic graph of x^2-29/16 is eight quarter-integers plus this one 3-cycle (no rational 6-cycle);")
print("      it is the unique centred monic quadratic over Q with a nondegenerate AP-supported 3-cycle; determinant-one")
print("      lift B with B^3=-I.  THM-4146-*: trace-one SL_2 lift, integral hexagon on X^2+2XY+13Y^2=48, the signed")
print("      Pythagorean template forces the 3:4:5 class and 29=5^2+2^2.  Rational 3-cycles of x^2+c: Poonen 1998 (CITED).")
print("(ii) 'the inverse square horizon (Z,Z^2) cancels the spatial translation after three steps': SCOPE, no map found;")
print("      the dynamics is x -> x^2 - 29/16, a single quadratic map, and its period-3 orbit is the content of THM-4139.")
# =============================================================================
hdr("P8. SECOND PART: Collatz peaks versus squares of 4k+1 (shortcut map C(n)=n/2, (3n+1)/2)")


def orbit_shortcut(n):
    out = [n]
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        out.append(n)
    return out


def orbit_full(n):
    out = [n]
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        out.append(n)
    return out


for seed in (5, 7, 23):
    os_ = orbit_shortcut(seed)
    of_ = orbit_full(seed)
    print("    seed %2d: shortcut orbit %s  peak %d" % (seed, os_, max(os_)))
    print("             full orbit peak %d" % max(of_))
check([max(orbit_shortcut(s)) for s in (5, 7, 23)] == [8, 26, 80], "shortcut peaks 8,26,80")
check([max(orbit_full(s)) for s in (5, 7, 23)] == [16, 52, 160], "full-map peaks 16,52,160")
print("(i) REFUTED ('endpoints 1, 26, 80'): the peaks of 5,7,23 are 8,26,80 (shortcut) or 16,52,160 (full map);")
print("      1 is the sink of every orbit, not a peak.  Incidentally 8=3^2-1, 26=5^2+1, 80=9^2-1 with bases 3,5,9.")
# structural fact: peaks of odd seeds are 2 mod 6 under the shortcut map
print("(ii) PROVED: for an odd seed n>1 the shortcut peak p is even (an odd p>1 is followed by (3p+1)/2>p), p=(3m+1)/2")
print("      with m odd and 3m+1 = 0 mod 4, so m = 1 mod 4 and p = 2 mod 6.  Under the full map the peak is 3m+1 = 4 mod 6.")
LIM = 10 ** 5
peaks = {}
peaks_full = {}
for n in range(3, LIM + 1, 2):
    o = orbit_shortcut(n)
    peaks[n] = max(o)
    m = n
    pf = n
    while m != 1:
        m = m // 2 if m % 2 == 0 else 3 * m + 1
        if m > pf:
            pf = m
    peaks_full[n] = pf
check(len(peaks) == 49999, "49,999 odd seeds in [3,10^5]")
check(all(p % 6 == 2 for p in peaks.values()), "all shortcut peaks are 2 mod 6")
check(all(p % 6 == 4 for p in peaks_full.values()), "all full-map peaks are 4 mod 6")
print("      FINITE-EXACT: all 49,999 shortcut peaks (odd seeds 3..10^5) are 2 mod 6; all full-map peaks are 4 mod 6.")
maxratio = max(Fraction(p, n) for n, p in peaks.items())
argmax = max(peaks, key=lambda n: Fraction(peaks[n], n))
print("      max peak/seed among them: %.3f at seed %d (peak %d)" % (float(maxratio), argmax, peaks[argmax]))
# hostile: 2^k-1
for k in (5, 10, 20):
    n = 2 ** k - 1
    o = orbit_shortcut(n)
    check(o[k] == 3 ** k - 1, "C^k(2^k-1) = 3^k-1")
print("      PROVED (hostile to 'growth bound by area'): C^j(2^k-1) = 3^j 2^{k-j} - 1 for j<=k, so the seed 2^k-1")
print("      climbs to 3^k-1 ~ n^{log2 3}; peak/seed is unbounded (checked k=5,10,20).")

# the statistic
def sq_window(p):
    """(name, base) if p = base^2 -1 or base^2 + 1 (p even so p != base^2); else None."""
    for q, name in ((p + 1, "s^2-1"), (p - 1, "s^2+1")):
        s_ = math.isqrt(q)
        if s_ * s_ == q:
            return name, s_
    return None


def v3(s_):
    k = 0
    while s_ % 3 == 0:
        s_ //= 3
        k += 1
    return k


hits = sum(1 for p in peaks.values() if (sq_window(p) is not None and sq_window(p)[1] % 4 == 1))
print("(iii) shortcut peaks in {(4k+1)^2-1,(4k+1)^2,(4k+1)^2+1}: %d of %d odd seeds" % (hits, len(peaks)))
check(hits == 248, "session lead's count 248 must be reproduced")
hits_any = sum(1 for p in peaks.values() if sq_window(p) is not None)
print("      peaks within +-1 of ANY square: %d (the window (4k+1)^2 itself is odd, never a peak)" % hits_any)
for r in range(4):
    hr = sum(1 for p in peaks.values() if sq_window(p) is not None and sq_window(p)[1] % 4 == r)
    print("      peaks within +-1 of a square of base = %d mod 4: %d" % (r, hr))
hits_full = sum(1 for p in peaks_full.values() if (sq_window(p) is not None and sq_window(p)[1] % 4 == 1))
print("      full-map peaks in the (4k+1)-square window: %d" % hits_full)
# seeds share peaks: count DISTINCT peak values and split hits into exact cells
from collections import Counter
cnt = Counter(peaks.values())
print("      distinct shortcut peak values: %d (top multiplicities %s)" % (len(cnt), cnt.most_common(4)))
check(len(cnt) == 20341, "20,341 distinct peaks")
cells = Counter()
cell_seeds = Counter()
for p, c in cnt.items():
    w = sq_window(p)
    if w is not None:
        key = (w[0], w[1] % 4, "3|s" if v3(w[1]) >= 1 else "3∤s")
        cells[key] += 1
        cell_seeds[key] += c
print("      PROVED: a peak p = 2 mod 6 with p = s^2-1 has s^2 = 3 mod 6 so 3|s; with p = s^2+1, s^2 = 1 mod 6 so 3∤s.")


def window_count(lo, hi, pred):
    c = 0
    for s_ in range(max(1, math.isqrt(max(lo - 2, 0))), math.isqrt(hi + 1) + 2):
        for q, name in ((s_ * s_ - 1, "s^2-1"), (s_ * s_ + 1, "s^2+1")):
            if lo <= q < hi and q % 6 == 2 and pred(s_, name):
                c += 1
    return c


dyadic = Counter(p.bit_length() - 1 for p in cnt)
null = {}
for name in ("s^2-1", "s^2+1"):
    for r in (1, 3):
        e = Fraction(0)
        for j, c in dyadic.items():
            lo, hi = 2 ** j, 2 ** (j + 1)
            class_size = sum(1 for m in range(lo + ((2 - lo) % 6), hi, 6))
            e += Fraction(c * window_count(lo, hi, lambda s_, nm, name=name, r=r: nm == name and s_ % 4 == r), class_size)
        null[(name, r)] = e
print("      cell table over DISTINCT peaks (null: a distinct peak is uniform on the class 2 mod 6 of its dyadic range):")
print("        window   base mod 4   3|s?   distinct peaks   seeds   null expected (distinct)")
for name in ("s^2-1", "s^2+1"):
    for r in (1, 3):
        d3 = "3|s" if name == "s^2-1" else "3∤s"
        key = (name, r, d3)
        print("        %-7s  %d            %-4s  %5d            %5d   %.1f" % (name, r, d3, cells[key], cell_seeds[key], float(null[(name, r)])))
check(cells[("s^2-1", 1, "3|s")] == 46 and cells[("s^2-1", 3, "3|s")] == 47, "s^2-1 cells 46 / 47")
check(cells[("s^2+1", 1, "3∤s")] == 49 and cells[("s^2+1", 3, "3∤s")] == 50, "s^2+1 cells 49 / 50")
check(cell_seeds[("s^2-1", 1, "3|s")] + cell_seeds[("s^2+1", 1, "3∤s")] == 248, "248 = 152 + 96")
plus_obs = cells[("s^2+1", 1, "3∤s")] + cells[("s^2+1", 3, "3∤s")]
plus_null = float(null[("s^2+1", 1)] + null[("s^2+1", 3)])
z_plus = (plus_obs - plus_null) / math.sqrt(plus_null)
minus_obs = cells[("s^2-1", 1, "3|s")] + cells[("s^2-1", 3, "3|s")]
minus_null = float(null[("s^2-1", 1)] + null[("s^2-1", 3)])
print("      +1 window: %d distinct vs null %.1f (Poisson z = %.2f): chance level." % (plus_obs, plus_null, z_plus))
print("      -1 window: %d distinct vs null %.1f: a real excess, split 46/47 between bases 1 and 3 mod 4." % (minus_obs, minus_null))
check(abs(z_plus) < 3, "+1 window at chance level")
check(minus_obs > 2 * minus_null, "-1 window excess is real")
# the mechanism: ladder a*2^k-1 -> a*3^k-1
explained = 0
unexplained = []
for p, c in cnt.items():
    w = sq_window(p)
    if w is None or w[0] != "s^2-1":
        continue
    found = False
    for k in range(1, 2 * v3(w[1]) + 1):
        if (p + 1) % 3 ** k == 0:
            a = (p + 1) // 3 ** k
            seed = a * 2 ** k - 1
            o = orbit_shortcut(seed)
            if seed % 2 == 1 and len(o) > k and o[k] == p and max(o) == p:
                found = True
                break
    if found:
        explained += 1
    else:
        unexplained.append(p)
print("      PROVED (ladder): C^j(a*2^k-1) = a*3^j*2^(k-j)-1 for j<=k; if a=m^2 and k is even the value after k odd")
print("        steps is (m*3^(k/2))^2-1, a square minus one.  All %d distinct '-1' peaks are of this form with the ladder" % minus_obs)
print("        seed a*2^k-1 (some k<=2v_3(s)) having that value as its peak: explained %d, unexplained %s." % (explained, unexplained))
check(explained == minus_obs and unexplained == [], "every -1 peak comes from the ladder")
print("      examples: 8=3^2-1 (seed 3=2^2-1; also 5=3*2-1), 80=9^2-1=3^4-1 (seed 15; also 53=27*2-1 on 23's orbit),")
print("        6560=3^8-1 (21 seeds), 59048=3^10-1 (48 seeds), 164024=25*3^8-1=405^2-1 (17 seeds), 4782968=3^14-1 (17 seeds).")
check(cnt[6560] == 21 and cnt[59048] == 48 and cnt[164024] == 17 and cnt[4782968] == 17, "ladder multiplicities")
print("      So the seed-weighted count 248 is inflated by shared ladder peaks; the 4k+1 selection is REFUTED (bases")
print("        1 and 3 mod 4 equally represented), and the only structure present is the ladder, which is exactly the")
print("        hostile of (ii) against 'growth bound by area'.")
# shifted windows (arbitrary residue control): {(4k+1)^2+s-1, +s, +s+1}
print("      shifted-window control, hits for peak in {(4k+1)^2+s-1,(4k+1)^2+s,(4k+1)^2+s+1}:")
row = []
for s in range(-12, 13, 2):
    hs = 0
    for p in peaks.values():
        for q in (p - s - 1, p - s, p - s + 1):
            if q >= 1:
                t = math.isqrt(q)
                if t * t == q and t % 4 == 1:
                    hs += 1
                    break
    row.append((s, hs))
print("        %s" % row)
print("      REFUTED ('accelerated Syracuse paths track the boundaries of perfect squares'): the hit count is at")
print("      chance level and shifted windows do as well; nothing selects the 4k+1 bases.")
# 196 and 169/160
primes12 = list(sympy.primerange(2, 40))[:12]
print("(iv) 196 = 14^2: sum of the first 12 primes %s = %d" % (primes12, sum(primes12)))
check(sum(primes12) == 197 and 196 == 14 ** 2, "sum of first 12 primes is 197")
print("      REFUTED ('196 = sum of the first 12 primes'): the sum is 197.  196 is not a peak (196 = 4 mod 6);")
o196 = orbit_shortcut(196)
print("      shortcut orbit of 196: %s (196 lies %d steps above 7; it is a predecessor, not a mirror)." % (o196[:o196.index(7) + 1], o196.index(7)))
check(o196[o196.index(7)] == 7 and o196.index(7) == 8, "196 reaches 7 in 8 shortcut steps")
seeds170 = sorted(n for n, p in peaks.items() if p == 170)
seeds168 = sorted(n for n, p in peaks.items() if p == 168)
seeds160 = sorted(n for n, p in peaks.items() if p == 160)
seeds160f = sorted(n for n, p in peaks_full.items() if p == 160)
print("      169 = 13^2: odd seeds <=10^5 with shortcut peak 170: %s; with peak 168: %s; with peak 160: %s" % (seeds170, seeds168, seeds160))
print("      full-map peak 160: seeds %s.  160 = 2^5*5 = 5*32 is on the trunk of 5 (160/2^5 = 5): TRUE arithmetic," % seeds160f)
print("      but 160 is the FULL-map peak of 23 (53 -> 160), its shortcut image being 80; 160 = 4 mod 6 can never be a")
print("      shortcut peak, and nothing relates 160 to 169 ('compression deficit' undefined): SCOPE / REFUTED.")
check(seeds160 == [] and 23 in seeds160f, "160 is a full-map peak of 23 and never a shortcut peak")
check(160 == 2 ** 5 * 5, "160 = 2^5 * 5")

# =============================================================================
hdr("P9. Descent certificate 2^K T^L(n) = 3^L n + B_L (inherited identity) and the open obligation")


def syracuse_word(n, L):
    """accelerated Syracuse: n odd -> (3n+1)/2^k; returns (k_1..k_L, final)."""
    ks = []
    for _ in range(L):
        m = 3 * n + 1
        k = 0
        while m % 2 == 0:
            m //= 2
            k += 1
        ks.append(k)
        n = m
    return ks, n


print("      n   L   k-word           K    3^L n + B      2^K T^L(n)   margin 2^K n-(3^L n+B)  descends")
for n, L in [(27, 5), (27, 41), (7, 2), (23, 3), (1, 1), (31, 5), (2 ** 20 - 1, 20)]:
    ks, fin = syracuse_word(n, L)
    K = sum(ks)
    Bv = sum(3 ** (L - 1 - j) * 2 ** sum(ks[:j]) for j in range(L))
    lhs = 2 ** K * fin
    rhs = 3 ** L * n + Bv
    check(lhs == rhs, "descent identity n=%d L=%d" % (n, L))
    margin = 2 ** K * n - rhs
    kw = str(ks) if len(ks) <= 8 else str(ks[:6])[:-1] + ", ...]"
    print("  %8d %3d  %-16s %4d  %-14s %-14s %-24s %s" % (n, L, kw, K, str(rhs)[:14], str(lhs)[:14], str(margin)[:24], margin > 0))
print("    PROVED (inherited, 05-knowledge/results/collatz_blueprint_20260921_synthesis.md section 6): with K_j the partial")
print("      sums of the halving word and B = sum_j 3^{L-1-j} 2^{K_j}, 2^K T^L(n) = 3^L n + B, hence")
print("      T^L(n) < n  iff  3^L n + B < 2^K n.  The paste's inequality is exactly this certificate.")
print("    The paste supplies no argument that every fixed n>1 eventually acquires a word with positive margin; 'peaks")
print("      run out of volume' is not a statement about any n.  The obligation remains OPEN, verbatim as in the audit.")
print("    27 reaches 1 after L=41 odd steps (K=%d, total 111 full-map steps)." % sum(syracuse_word(27, 41)[0]))
# first L with descent for 27
ks_all, _ = syracuse_word(27, 60)
firstL = None
for L in range(1, 61):
    K = sum(ks_all[:L])
    Bv = sum(3 ** (L - 1 - j) * 2 ** sum(ks_all[:j]) for j in range(L))
    if 2 ** K * 27 - (3 ** L * 27 + Bv) > 0:
        firstL = L
        break
K37 = sum(ks_all[:firstL])
check(firstL == 37 and K37 == 59 and syracuse_word(27, 37)[1] == 23, "27 first descends below itself at L=37, K=59, value 23")
print("      FINITE-EXACT: for n=27 the first L with positive margin is L=%d (K=%d, T^37(27)=23; 37+59=96 is the" % (firstL, K37))
print("        classical total stopping time of 27 under the full map).")

# =============================================================================
hdr("P10. Claim 6 and the remaining phrases")
wheel = [r for r in range(30) if math.gcd(r, 30) == 1]
check(wheel == [1, 7, 11, 13, 17, 19, 23, 29] and len(wheel) == 8, "mod-30 wheel")
print("    mod-30 wheel: reduced classes %s (phi(30)=8) -- the only defined object in 'mod 30 primes framework'." % wheel)
print("    'Unified Arithmetic Cycle {f,+,*,g}', 'modular stabilizer q = x (mod AB)', 'macro-growth fractures back into")
print("      local modular grids', 'Bott periodicity transforming tetration into pentation', 'quasicrystalline Fourier")
print("      transform matrix': SCOPE -- no definitions, no map found.  (Bott periodicity is the 8-fold periodicity of")
print("      the stable homotopy of O and U (CITED: Bott 1959); nothing connects it to hyperoperations.)")

# =============================================================================
hdr("SUMMARY TABLE")
rows = [
    ("1 core identity T((A+1)(B+1)-1)-T(AB-1)=T(A)+T(B)+Shear", "FALSE / TRUE-with-repair", "A=B=1: 6 vs 2; repair: shift both arguments by -1 (P1 iii)"),
    ("1 horizontal/vertical flows", "TRUE-with-repair", "S=AB(A+B)+C(A,2)+C(B,2)=paste's RHS, not its LHS (P1 iv)"),
    ("1 'shear = swept volume', 'conservation law'", "SCOPE", "true shear AB(A+B+2)=2A T(B)+2B T(A); no volume defined"),
    ("2 tetrahedral law", "TRUE", "Vandermonde split (2,1) (P2)"),
    ("2 pentatope law", "TRUE", "Vandermonde split (2,2) (P2)"),
    ("2 'prisms locking against a plate'", "SCOPE", "a set-partition count; no shear"),
    ("3 AgB=AB*(AgB)", "FALSE", "forces AgB=0 unless AB=1"),
    ("3 AgB=(AB)^(AB) and 1gB=B!", "FALSE", "B=2: 4 vs 2; consistent reading A!B! (P3)"),
    ("3 AgA=(A!)^2", "TRUE-with-repair", "holds for AgB=A!B! only"),
    ("3 Pell / square-triangular mapping", "SCOPE", "no map found"),
    ("4 free edges = T(n-2)", "TRUE", "C(n,2)-(n-1)=C(n-1,2)"),
    ("4 path edges carry 0 bits", "FALSE", "code costs C(n,2)+log2(n!/2^(n-1)) bits; sum h = n! 2^T(n-2) (P4)"),
    ("4 Wythoff diffraction / genus", "SCOPE", "no object"),
    ("5 horizon lines", "SCOPE", "definitions of diagonals"),
    ("5 63 bottleneck", "TRUE-with-repair", "63=2^6-1 is Bang's unique exception via 9=3^2 (P5)"),
    ("5 341 bottleneck", "FALSE", "4^5-1 has primitive primes 11,31; Cipolla pseudoprime on the R-trunk (P6)"),
    ("5 rational 3-cycle", "TRUE (content = THM-4139/4146)", "'horizon cancels translation': SCOPE (P7)"),
    ("6 unified cycle, Bott, quasicrystal FT", "SCOPE", "mod-30 wheel only (P10)"),
    ("II peaks 1,26,80", "FALSE", "peaks 8,26,80 (shortcut) / 16,52,160 (full) (P8)"),
    ("II peaks track (4k+1)^2", "FALSE", "bases 1/3 mod 4 equally hit (46/47, 49/50 distinct); -1 excess = ladder a*3^k-1 (P8)"),
    ("II 196 = sum of first 12 primes", "FALSE", "sum = 197"),
    ("II 160 = 2^5*5 on trunk of 5", "TRUE", "160 is 23's full-map peak; unrelated to 169"),
    ("II 3^L n + B_L < 2^K n", "TRUE (inherited)", "descent certificate; obligation OPEN (P9)"),
    ("II 'peaks bound by area of squares' proves descent", "FALSE", "2^k-1 climbs to 3^k-1; no proof (P8,P9)"),
]
for c, v, w in rows:
    print("  | %-58s | %-30s | %s" % (c, v, w))

print()
print("FAILS: %d" % len(FAILS))
print("ALL CHECKS PASSED" if not FAILS else "SOME CHECKS FAILED")
