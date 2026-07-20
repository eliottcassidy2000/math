"""opus-2026-07-20-S420: SHEAR CALCULUS on the n*2^x+1 grid + the census of
triangular-number continuations.

Owner: shear the (x,n) grid of n*2^x+1 into a triangle (columns down 1) or a
Fibonacci-analogue (down 2); sum or MULTIPLY the sheared rows; find as many
meaningful continuations of the triangular numbers as possible; compare all
under shears 0/1/2(/3); triangular numbers are relational (edges of K_n).

Machinery: shear_s row m of grid A = { A(x, m - s*x) : x >= 0, m - s*x >= 0 }.
Exact minimal-recurrence detector (Fraction Gaussian elimination, no floats).
"""
from fractions import Fraction
from math import comb, gcd, isqrt

# ---------- exact minimal linear recurrence detector ----------
def minrec(seq, maxord=8):
    L = len(seq)
    for r in range(1, maxord + 1):
        if L < 2 * r + 3:
            break
        # solve first solvable r x r window exactly, then verify globally
        for start in range(r, L - r):
            M = [[Fraction(seq[m - i]) for i in range(1, r + 1)] + [Fraction(seq[m])]
                 for m in range(start, start + r)]
            # gaussian elimination
            ok = True
            for col in range(r):
                piv = next((i for i in range(col, r) if M[i][col] != 0), None)
                if piv is None:
                    ok = False
                    break
                M[col], M[piv] = M[piv], M[col]
                M[col] = [v / M[col][col] for v in M[col]]
                for i in range(r):
                    if i != col and M[i][col] != 0:
                        M[i] = [a - M[i][col] * b for a, b in zip(M[i], M[col])]
            if not ok:
                continue
            cs = [M[i][r] for i in range(r)]
            if all(sum(cs[i] * seq[m - i - 1] for i in range(r)) == seq[m]
                   for m in range(r, L)):
                return r, cs
    return None, None

def rec_str(r, cs):
    if r is None:
        return "NONE (no constant-coefficient recurrence up to order 8 -- honest)"
    return "a(m) = " + " + ".join(f"({c})a(m-{i+1})" for i, c in enumerate(cs) if c != 0)

# ---------- grids: continuations of the triangular numbers ----------
def qbin2(n, k):
    if k < 0 or k > n:
        return 0
    num = den = 1
    for i in range(k):
        num *= 2 ** (n - i) - 1
        den *= 2 ** (i + 1) - 1
    return num // den

def stirling2(n, k):
    if k < 0 or k > n:
        return 0
    if n == 0:
        return 1 if k == 0 else 0
    return k * stirling2(n - 1, k) + stirling2(n - 1, k - 1)

GRIDS = {
    "proth n*2^x+1":      lambda x, n: n * 2 ** x + 1,
    "pure n*2^x":         lambda x, n: n * 2 ** x,
    "pascal C(n,x)":      lambda x, n: comb(n, x) if x <= n else 0,
    "simplicial C(n+x,x+1)": lambda x, n: comb(n + x, x + 1),
    "faulhaber sum j^x":  lambda x, n: sum(j ** x for j in range(1, n + 1)),
    "polygonal P(x+3,n)": lambda x, n: ((x + 1) * n * n - (x - 1) * n) // 2,
    "centered (x+3)T(n-1)+1": lambda x, n: (x + 3) * (n - 1) * n // 2 + 1 if n >= 1 else 1,
    "qgauss2 [n,x]_2":    lambda x, n: qbin2(n, x),
    "stirling2 S(n,x)":   lambda x, n: stirling2(n, x),
}

def shear_row(A, s, m, square_window=None):
    if s == 0:
        return [A(x, m) for x in range(0, (square_window or m) + 1)]
    return [A(x, m - s * x) for x in range(0, m // s + 1)]

print("=" * 88)
print("(1) SHEAR SPECTRUM: row sums of every grid at shears s = 0(window x<=m), 1, 2, 3")
print("=" * 88)
NAMED = []
for name, A in GRIDS.items():
    for s in (0, 1, 2, 3):
        MMAX = 22 if name != "qgauss2 [n,x]_2" else 13
        seq = [sum(shear_row(A, s, m, square_window=m)) for m in range(0, MMAX)]
        r, cs = minrec(seq)
        tag = ""
        if seq[1:9] == [2 ** (m + 1) - 1 for m in range(1, 9)]:
            tag = "  == MERSENNE 2^(m+1)-1 !!"
        if seq[1:9] == [2 ** m - 1 for m in range(1, 9)]:
            tag = "  == MERSENNE 2^m-1 !!"
        if seq[2:9] == [2 ** (m - 1) for m in range(2, 9)]:
            tag = "  == powers of 2"
        fib = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
        if any(seq[2:9] == fib[j:j + 7] for j in range(4)):
            tag = "  == FIBONACCI !!"
        if seq[2:10] == [fib[m + 1] - 1 for m in range(1, 9)]:
            tag = "  == Fibonacci - 1"
        print(f"  {name:26s} s={s}: {seq[:10]}")
        print(f"      {'':22s} rec: {rec_str(r, cs)}{tag}")
        if tag:
            NAMED.append((name, s, tag.strip()))
print("  NAMED HITS:", NAMED)

print()
print("=" * 88)
print("(2) PROTH TRIANGLE (shear-1) HEADLINES")
print("=" * 88)
P = GRIDS["proth n*2^x+1"]
print("  rows m=0..5:", [shear_row(P, 1, m) for m in range(6)])
print("  row sums == 2^(m+1)-1 for m<=20:",
      all(sum(shear_row(P, 1, m)) == 2 ** (m + 1) - 1 for m in range(21)),
      "(proof: sum n*2^x over x+n=m is A000295(m+1) = 2^(m+1)-m-2; +(m+1) ones)")
# collision law: n*2^x = n'*2^x' within a row <=> n(2^d - 1) = d
sols = [(n, d) for d in range(1, 60) for n in range(1, 200)
        if n * (2 ** d - 1) == d]
print("  collision law n(2^d-1)=d, all solutions n,d>=1:", sols,
      "-> UNIQUE (n=1,d=1): every row m>=2 holds the FERMAT number 2^(m-1)+1 exactly twice")
prods = []
for m in range(0, 13):
    row = shear_row(P, 1, m)
    pr = 1
    for v in row:
        pr *= v
    sq = isqrt(pr) ** 2 == pr
    prods.append((m, pr if m <= 6 else f"~2^{pr.bit_length()-1}", sq))
print("  row products (m, value, perfect square?):", prods)
print("  -> squares ONLY at m=2 (9=3^2) and m=3 (100=10^2); Fermat square (2^(m-1)+1)^2 divides every row product")
import math
scale = [(m, round(math.log2(math.prod(shear_row(P, 1, m))) / comb(m, 2), 4))
         for m in (6, 10, 14, 18, 22)]
print("  log2(row product)/C(m,2) ->", scale, "-> 1: PRODUCTS live on the TOURNAMENT scale 2^C(m,2);")
print("     sums live on the subset scale 2^m.  Sum counts elements; product counts orientations.")
pp = []
for m in range(2, 15):
    row = shear_row(P, 1, m)
    cnt = 0
    for v in row:
        if v > 1 and all(v % p for p in range(2, isqrt(v) + 1)):
            cnt += 1
    pp.append(cnt)
print("  primes per shear-1 row (m=2..14) [Proth-prime census]:", pp)

print()
print("=" * 88)
print("(3) SHEAR-2 = THE FIBONACCI MOVE: the growth dial")
print("=" * 88)
s2 = [sum(shear_row(P, 2, m)) for m in range(24)]
r, cs = minrec(s2)
print("  proth shear-2 sums:", s2[:12])
print("  minimal recurrence:", rec_str(r, cs))
print("  closed form: even m=2t: 2^(t+2)-t-3; odd m=2t+1: 3*2^(t+1)-t-4; growth (sqrt2)^m:",
      all((s2[m] == 2 ** (m // 2 + 2) - m // 2 - 3 if m % 2 == 0
           else s2[m] == 3 * 2 ** ((m - 1) // 2 + 1) - (m - 1) // 2 - 4) for m in range(24)))
print("  char poly (x^2-2)(x-1)^2(x+1): roots +-sqrt2 (two-phase), 1 (double), -1")
print("  GROWTH DIAL: Pascal shear-s -> Pisot root of x^s = x^(s-1)+1 (2, phi, Narayana-cows...);")
print("               Proth  shear-s -> 2^(1/s) exactly (2, sqrt2, cbrt2...): exponential grids dial RADICALS,")
print("               combinatorial grids dial PISOT numbers.  Both -> 1 as s -> infinity.")
print("  charm: proth shear-2 sum at m=11 =", s2[11], "(= Phi_6(14) = the deep well denominator, charm-grade)")

print()
print("=" * 88)
print("(4) THE MUSEUM OF IMPERSONATIONS (low-order jets colliding -- the session theme)")
print("=" * 88)
lazy = [1 + m * (m + 1) // 2 for m in range(10)]
print("  A. proth shear-2 sums vs lazy caterer 1+T(m):", s2[:8], "vs", lazy[:8],
      "-> split at m=5 (18 vs 16)")
moser = [comb(n - 1, 4) + comb(n - 1, 2) + 1 for n in range(1, 11)]
print("  B. Moser circle vs 2^(n-1):", moser[:8], "-> split at n=6 (31 vs 32)")
cp = [1, 3]
while len(cp) < 8:
    cp.append(2 * cp[-1] + cp[-2])
proth_diag = [k * 2 ** (k - 2) + 1 for k in range(2, 9)]
print("  C. Proth diagonal k*2^(k-2)+1 vs companion Pell:", proth_diag[:6], "vs", cp[1:7],
      "-> split at 5th term (97 vs 99)")
F = GRIDS["faulhaber sum j^x"]
diag = [sum(shear_row(F, 2, m)) for m in range(1, 10)]
pred = [None] * len(diag)
brk = None
for i in range(4, len(diag)):
    p = diag[i - 1] + diag[i - 2] + diag[i - 4]
    if p != diag[i]:
        brk = (i + 1, p, diag[i])
        break
print("  D. MY OWN S419 MIRAGE: faulhaber shear-2 sums", diag,
      f"; law a(n)=a(n-1)+a(n-2)+a(n-4) BREAKS at m={brk[0]} (predicts {brk[1]}, actual {brk[2]})")
print("     -> S419's 'verified on all 3 instances' was a category-E jet illusion. MISTAKE-191 filed.")
print("  E. repo's own: width of G_n = C(n-2,floor) -- exact n=3..6, fails n=7 (predicted 10, actual 15).")
print("  MECHANISM: all five are partial low-order jets of different growth classes (quadratic /")
print("  exponential / Pisot / radical); the shear calculus names the true class and predicts the split.")

print()
print("=" * 88)
print("(5) THE RELATIONAL SYNTHESIS: T_n = n in the relational sense")
print("=" * 88)
print("  2*T_m + 1 = m^2+m+1 = Phi_6(m+1) = |PG(2,m)|:",
      all(2 * (m * (m + 1) // 2) + 1 == m * m + m + 1 for m in range(1, 200)))
print("  deep well: 183 = 2*T_13 + 1 =", 2 * (13 * 14 // 2) + 1,
      "= |PG(2,13)| -- the LRC deep well denominator IS the doubled-relation-plus-one of 13.")
print("  the owner's (1,n)=2n+1 column, evaluated AT triangular arguments, = projective plane sizes.")
print("  tournament tower: #labeled tournaments on n = 2^T(n-1); staircase tile count m = C(n-1,2) = T(n-2):")
print("    the repo's tiling hypercube Q_m has DIMENSION a triangular number -- relation count as exponent.")
print("  centered polygonals = s*T(n-1) + 1: the bilinear m*b+1 rim WITH TRIANGULAR BASE.")
print("  hexagonal fold-back: H(n) = T(2n-1):",
      all(n * (2 * n - 1) == (2 * n - 1) * 2 * n // 2 for n in range(1, 50)),
      "(the polygonal tower folds into the triangular spine)")
print()
print("  CONTINUATION AXES CENSUS (the wide berth): 8 independent axes out of T_n:")
print("   1. arity      : C(n,2) -> C(n,k)          [k-ary relations, hypergraph edges]")
print("   2. geometry   : P(3,n) -> P(s,n)          [polygonal; hexagonal folds back to T]")
print("   3. dimension  : C(n+1,2) -> C(n+d-1,d)    [simplicial: tetrahedral, pentatope...]")
print("   4. exponent   : sum j^1 -> sum j^k        [Faulhaber tower, S419's third perspective]")
print("   5. field      : C(n,2) -> [n,2]_q         [q-analog: 2-subspaces of F_q^n; q=2 the dyadic shadow]")
print("   6. orientation: T_n -> 2^(T_n)            [tournaments; sums=elements, products=orientations]")
print("   7. centering  : T_n -> s*T(n-1)+1         [centered polygonal = bilinear rim on base T]")
print("   8. projective : T_n -> 2T_n+1 = |PG(2,n)| [doubling+duty: the deep well axis]")
print("  q-check [n,2]_2 (dyadic triangulars):", [qbin2(n, 2) for n in range(2, 8)],
      "(1,7,35,155,651,2667: 2-subspaces of F_2^n)")
