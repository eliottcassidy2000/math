#!/usr/bin/env python3
"""
triangular_continuations_wide_birth_boxeph_S148.py  (HYP-8175)

THE WIDE BIRTH: every meaningful continuation of the triangular numbers we can
build, sheared down 0/1/2/3, summed / alternated / multiplied, and every derived
sequence pushed through an exact recognizer (named-sequence library + minimal
linear recurrence + polynomiality).  Owner directive S148.

ANCHOR (owner): "the triangular numbers arise from relation itself, the edges of
complete graphs; the nth triangular number is n in the relational sense."
T(n) = C(n,2) = #edges(K_n) = #relations among n things.  The continuations are
the ways to deform "count the pair-structures":

  FAMILY        G(x, n)                 triangular section     axis meaning
  A pascal      C(n+x, x)               x = 2 diag             dimension (simplicial)
  B polygonal   P_{3+x}(n)              x = 0                  step (gnomon size)
  C pyramidal   Y_{3+x}(n)              --  (contains sq-pyr)  step, then stack
  D faulhaber   S_x(n) = sum k^x        x = 1                  exponent (S147 3rd view)
  E qbinom      [n choose 2]_{q=1+x}    x = 0                  ground field (q-deform)
  F proth       n*2^x + 1               --  (owner grid)       tower height
  G subsets     C(n, x)                 x = 2 col              #groups, subset sense
  H partitions  S(n, x)  (Stirling 2)   x = 2 col = 2^{n-1}-1  #groups, partition sense
  I cycles      c(n, x)  (Stirling 1)   x = 2 col = (n-1)!H    #groups, cycle sense
  J owner       T(m,j) via S147 form    col 2                  the riddle triangle
  K moser       R_x(n) = sum_{i<=x}C(n,i)  x = 2 (lazy caterer-1?)  truncation (Moser/2^n break)

THE RELATIONAL TRIAD (G,H,I at x = 2): the three fundamental "organize n things
into TWO groups" counts -- 2-subsets C(n,2) = TRIANGULAR, 2-blocks S(n,2) =
MERSENNE 2^{n-1}-1, 2-cycles c(n,2) = (n-1)! * H_{n-1} = HARMONIC numerators.
All three surfaced independently in the owner's S147 riddle (triangular column,
Mersenne-Moser break 1,3,7,15,29, harmonic series).  Verified exactly below.

SHEAR-t: column x moved down t*x rows: sheared row m = { G(x, m - t*x) }.
t=1 = Pascal-style triangle; t=2 = Fibonacci-style (shallow diagonals); t=0 =
the square itself (finite rows only).  CONTROL: pascal shear-1 sums = 2^m,
shear-2 sums = Fibonacci (asserted).  NEW EXACT LAWS (derived by hand, asserted):
  * proth shear-1 row sums = 2^{m+1} - 1        (MERSENNE -- the S147 break!)
  * proth shear-2 row sums: m=2K   -> 2^{K+2} - K - 3
                            m=2K+1 -> 3*2^{K+1} - K - 4   (Lucas-flavored start 1,2,4,7,11,18)

boxeph-2026-07-20-S148.  Pure Python, exact integers/Fractions.
"""

from fractions import Fraction as Fr
from math import comb, prod, gcd
from functools import lru_cache

M_MAX = 16          # sheared rows computed for m = 0..M_MAX
SEQ_LEN = 15        # terms fed to the recognizer
PROD_MAX = 8        # row products only up to this m (they explode)

# ----------------------------------------------------------------------------
# the families: g(x, n) -> int or None (out of domain).  n is the ROW index the
# shear consumes; x is the COLUMN being moved down.
# ----------------------------------------------------------------------------

@lru_cache(maxsize=None)
def stir2(n, k):
    if n == 0: return 1 if k == 0 else 0
    if k == 0: return 0
    return k * stir2(n-1, k) + stir2(n-1, k-1)

@lru_cache(maxsize=None)
def stir1(n, k):   # unsigned cycle numbers
    if n == 0: return 1 if k == 0 else 0
    if k == 0: return 0
    return (n-1) * stir1(n-1, k) + stir1(n-1, k-1)

def fam_pascal(x, n):     return comb(n + x, x) if n >= 0 else None
def fam_polygonal(x, n):  # k-gonal, k = 3+x, n >= 0
    if n < 0: return None
    k = 3 + x
    return ((k - 2) * n * n - (k - 4) * n) // 2
def fam_pyramidal(x, n):  # k-gonal pyramid, k = 3+x
    if n < 0: return None
    k = 3 + x
    return sum(((k - 2) * j * j - (k - 4) * j) // 2 for j in range(n + 1))
def fam_faulhaber(x, n):  # S_x(n) = 1^x + ... + n^x ; S_0(n) = n
    if n < 0: return None
    return sum(j ** x for j in range(1, n + 1))
def fam_qbinom(x, n):     # [n choose 2]_q at q = 1+x  (q=1: triangular C(n,2))
    if n < 0: return None
    q = 1 + x
    if n < 2: return 0
    if q == 1: return comb(n, 2)
    num = (q ** n - 1) * (q ** (n - 1) - 1)
    den = (q * q - 1) * (q - 1)
    assert num % den == 0
    return num // den
def fam_proth(x, n):      return n * 2 ** x + 1 if n >= 0 else None
def fam_subsets(x, n):    return comb(n, x) if 0 <= x <= n else None
def fam_partitions(x, n): return stir2(n, x) if 0 <= x <= n else None
def fam_cycles(x, n):     return stir1(n, x) if 0 <= x <= n else None
def _S(p, n): return sum(k ** p for k in range(1, n + 1))
def fam_owner(x, n):      # owner triangle T(m,j): m = n, j = x+1 (columns as x)
    m, j = n, x + 1
    if not (1 <= j <= m): return None
    c = comb(m - 6, j - 4) if 0 <= j - 4 <= m - 6 else 0
    return _S(j - 1, m - j + 1) + c
def fam_moser(x, n):      # partial Pascal row sums R_x(n) = sum_{i<=x} C(n,i)
    if n < 0: return None
    return sum(comb(n, i) for i in range(0, x + 1))

FAMS = [
    ("pascal",     fam_pascal,     False),
    ("polygonal",  fam_polygonal,  False),
    ("pyramidal",  fam_pyramidal,  False),
    ("faulhaber",  fam_faulhaber,  False),
    ("qbinom2",    fam_qbinom,     False),
    ("proth",      fam_proth,      False),
    ("subsets",    fam_subsets,    True),   # True = rows naturally finite -> shear-0 ok
    ("partitions", fam_partitions, True),
    ("cycles",     fam_cycles,     True),
    ("owner",      fam_owner,      True),
    ("moser",      fam_moser,      False),
]

# ----------------------------------------------------------------------------
# shear machinery
# ----------------------------------------------------------------------------

def sheared_row(g, t, m, xmax=64):
    """row m of the shear-t triangle: [G(x, m - t*x)] over the domain."""
    row = []
    for x in range(0, xmax + 1):
        n = m - t * x
        if t > 0 and n < 0: break
        v = g(x, n)
        if v is None:
            if t == 0: break
            continue
        row.append(v)
    return row

# ----------------------------------------------------------------------------
# recognizer: named library, minimal linear recurrence, polynomiality
# ----------------------------------------------------------------------------

def _mk(f, k=24): return [f(i) for i in range(k)]
def _fib():
    a, b = 0, 1; out = []
    for _ in range(26): out.append(a); a, b = b, a + b
    return out
def _lucas():
    a, b = 2, 1; out = []
    for _ in range(26): out.append(a); a, b = b, a + b
    return out
def _pell():
    a, b = 0, 1; out = []
    for _ in range(22): out.append(a); a, b = b, 2 * b + a
    return out
def _jacobsthal():
    a, b = 0, 1; out = []
    for _ in range(24): out.append(a); a, b = b, b + 2 * a
    return out
def _trib():
    s = [0, 0, 1]
    while len(s) < 24: s.append(s[-1] + s[-2] + s[-3])
    return s
def _padovan():
    s = [1, 1, 1]
    while len(s) < 26: s.append(s[-2] + s[-3])
    return s
def _catalan(): return [comb(2 * i, i) // (i + 1) for i in range(16)]
def _motzkin():
    m = [1, 1]
    for i in range(2, 18):
        m.append(((2 * i + 1) * m[-1] + 3 * (i - 1) * m[-2]) // (i + 2))
    return m
def _bell():
    b = [[1]]
    for i in range(1, 16):
        row = [b[-1][-1]]
        for v in b[-1]: row.append(row[-1] + v)
        b.append(row)
    return [r[0] for r in b]
def _a000522(): # sum n!/k!
    out = []
    from math import factorial
    for n in range(14):
        out.append(sum(factorial(n) // factorial(k) for k in range(n + 1)))
    return out

from math import factorial
LIB = {
    "Fibonacci": _fib(), "Lucas": _lucas(), "Pell": _pell(),
    "Jacobsthal": _jacobsthal(),
    "Jacobsthal-Lucas": _mk(lambda i: 2 ** i + (-1) ** i * 1 if i else 2, 20),
    "Mersenne 2^n-1": _mk(lambda i: 2 ** i - 1, 22),
    "2^n+1": _mk(lambda i: 2 ** i + 1, 20),
    "2^n": _mk(lambda i: 2 ** i, 22),
    "3^n": _mk(lambda i: 3 ** i, 16),
    "(3^n+1)/2": _mk(lambda i: (3 ** i + 1) // 2, 16),
    "(3^n-1)/2": _mk(lambda i: (3 ** i - 1) // 2, 16),
    "tribonacci": _trib(), "Padovan": _padovan(),
    "Catalan": _catalan(), "Motzkin": _motzkin(), "Bell": _bell(),
    "factorial": _mk(lambda i: factorial(i), 14),
    "triangular": _mk(lambda i: i * (i + 1) // 2, 22),
    "squares": _mk(lambda i: i * i, 22),
    "pronic": _mk(lambda i: i * (i + 1), 22),
    "tetrahedral": _mk(lambda i: comb(i + 2, 3), 20),
    "Cullen n*2^n+1": _mk(lambda i: i * 2 ** i + 1, 18),
    "Woodall n*2^n-1": _mk(lambda i: i * 2 ** i - 1, 18),
    "A000522 sum n!/k!": _a000522(),
    "lazy caterer": _mk(lambda i: i * (i + 1) // 2 + 1, 22),
    "Moser circle": _mk(lambda i: 1 + comb(i, 2) + comb(i, 4), 20),
    "Eulerian 2^n-n-1": _mk(lambda i: 2 ** i - i - 1, 20),
    "2^n-n": _mk(lambda i: 2 ** i - i, 20),
    "n*2^(n-1)": _mk(lambda i: i * 2 ** (i - 1) if i else 0, 18),
}

def lib_match(seq):
    """contiguous subsequence match against the library, offset 0..4."""
    hits = []
    s = [v for v in seq if v is not None]
    if len(s) < 5: return hits
    for name, ref in LIB.items():
        for off in range(0, 5):
            if off + len(s) <= len(ref) and ref[off:off + len(s)] == s:
                hits.append("%s[%d:]" % (name, off)); break
            # also allow matching a truncation (>= 7 terms)
            L = min(len(s), len(ref) - off)
            if L >= 7 and ref[off:off + L] == s[:L]:
                hits.append("%s[%d:]~%dt" % (name, off, L)); break
    return hits

def min_linrec(seq, maxord=6):
    """minimal exact constant-coeff linear recurrence; needs >= 2r+3 terms."""
    s = [Fr(v) for v in seq]
    for r in range(1, maxord + 1):
        if len(s) < 2 * r + 3: break
        # solve a(i) = sum_j c_j a(i-j), i = r..len-1 (exact Gauss)
        rows = [[s[i - j] for j in range(1, r + 1)] + [s[i]] for i in range(r, len(s))]
        nc = r; A = [row[:] for row in rows]; piv = []; rk = 0
        for c in range(nc):
            p = next((i for i in range(rk, len(A)) if A[i][c] != 0), None)
            if p is None: continue
            A[rk], A[p] = A[p], A[rk]
            pr = A[rk]
            for i in range(len(A)):
                if i != rk and A[i][c] != 0:
                    f = A[i][c] / pr[c]
                    for cc in range(nc + 1): A[i][cc] -= f * pr[cc]
            piv.append(c); rk += 1
        if any(all(row[c] == 0 for c in range(nc)) and row[nc] != 0 for row in A):
            continue
        sol = [Fr(0)] * nc
        for i, c in enumerate(piv): sol[c] = A[i][nc] / A[i][c]
        if all(s[i] == sum(sol[j - 1] * s[i - j] for j in range(1, r + 1))
               for i in range(r, len(s))):
            return r, sol
    return None, None

def poly_degree(seq, maxdeg=7):
    s = [Fr(v) for v in seq]
    for d in range(0, maxdeg + 1):
        t = s[:]
        for _ in range(d + 1):
            t = [t[i + 1] - t[i] for i in range(len(t) - 1)]
        if len(t) >= 2 and all(v == 0 for v in t):
            return d
    return None

def identify(seq):
    """one-line identification string for a sequence."""
    if len(set(seq)) == 1: return "constant %s" % seq[0]
    tags = lib_match(seq)
    pd = poly_degree(seq)
    if pd is not None: tags.append("POLY deg %d" % pd)
    else:
        r, coef = min_linrec(seq)
        if r is not None:
            tags.append("linrec r=%d %s" % (r, [str(c) for c in coef]))
    return "; ".join(tags) if tags else "UNRECOGNIZED"

# ----------------------------------------------------------------------------
# (0) hand-derived exact laws + controls, asserted
# ----------------------------------------------------------------------------
print("=" * 100)
print("(0) CONTROLS AND NEW EXACT LAWS (asserted)")
print("=" * 100)
ps1 = [sum(sheared_row(fam_pascal, 1, m)) for m in range(12)]
ps2 = [sum(sheared_row(fam_pascal, 2, m)) for m in range(14)]
assert ps1 == [2 ** m for m in range(12)], ps1
fib = _fib()
assert ps2 == fib[1:15], (ps2, fib[1:15])
print("CONTROL pascal: shear-1 sums = 2^m OK ; shear-2 sums = Fibonacci OK")
pr1 = [sum(sheared_row(fam_proth, 1, m)) for m in range(14)]
assert pr1 == [2 ** (m + 1) - 1 for m in range(14)], pr1
print("LAW 1 (NEW): proth shear-1 row sums = 2^(m+1)-1 = MERSENNE   %s..." % pr1[:8])
pr2 = [sum(sheared_row(fam_proth, 2, m)) for m in range(16)]
def proth_s2(m):
    K = m // 2
    return 2 ** (K + 2) - K - 3 if m % 2 == 0 else 3 * 2 ** (K + 1) - K - 4
assert pr2 == [proth_s2(m) for m in range(16)], pr2
print("LAW 2 (NEW): proth shear-2 sums: even m=2K: 2^(K+2)-K-3 ; odd m=2K+1: 3*2^(K+1)-K-4")
print("             = %s...  (Lucas 1,2,4,7,11,18 start, then splits: ~C*sqrt(2)^m)" % pr2[:10])
print("             recognizer says:", identify(pr2))
# the relational triad
H = [sum(Fr(1, k) for k in range(1, n + 1)) for n in range(0, 12)]
triad_ok = all(fam_subsets(2, n) == n * (n - 1) // 2 for n in range(2, 12)) \
    and all(fam_partitions(2, n) == 2 ** (n - 1) - 1 for n in range(2, 12)) \
    and all(Fr(fam_cycles(2, n + 1), factorial(n)) == H[n] for n in range(1, 11))
print("RELATIONAL TRIAD: C(n,2)=triangular, S(n,2)=2^(n-1)-1 MERSENNE, c(n+1,2)/n! = H_n :",
      "ALL VERIFIED" if triad_ok else "FAIL")
assert triad_ok
q2 = [fam_qbinom(1, n) for n in range(0, 10)]
print("q-TRIANGULAR at q=2 (2-subspaces of F2^n, the BINARY relational count):", q2)
print("  linrec:", identify(q2[2:]))

# ----------------------------------------------------------------------------
# (1) the full comparison: family x shear -> sums / altsums identified
# ----------------------------------------------------------------------------
print("\n" + "=" * 100)
print("(1) COMPARISON MATRIX: row sums of each shear (t=0 only for finite-row families)")
print("=" * 100)
matrix_notes = []
for name, g, finite in FAMS:
    shears = ([0] if finite else []) + [1, 2, 3]
    for t in shears:
        sums, alts = [], []
        for m in range(M_MAX + 1):
            row = sheared_row(g, t, m)
            if not row: continue
            sums.append(sum(row))
            alts.append(sum((-1) ** i * v for i, v in enumerate(row)))
        sums = sums[:SEQ_LEN]; alts = alts[:SEQ_LEN]
        idn = identify(sums)
        ida = identify(alts)
        print("%-11s t=%d  SUM %-22s -> %s" % (name, t, str(sums[:7]) + "...", idn))
        print("%-11s      ALT %-22s -> %s" % ("", str(alts[:7]) + "...", ida))
        matrix_notes.append((name, t, idn, ida))

# ----------------------------------------------------------------------------
# (2) products (owner: "or in fact take their product")
# ----------------------------------------------------------------------------
print("\n" + "=" * 100)
print("(2) ROW PRODUCTS (m <= %d): growth signatures + any small-sequence hits" % PROD_MAX)
print("=" * 100)
for name, g, finite in FAMS:
    for t in ([0] if finite else []) + [1, 2]:
        prods = []
        for m in range(PROD_MAX + 1):
            row = [v for v in sheared_row(g, t, m) if v != 0]
            if not row: continue
            prods.append(prod(row))
        if len(prods) >= 5:
            print("%-11s t=%d  PROD %s" % (name, t, prods[:7]))
            hits = lib_match(prods)
            if hits: print("%-11s      -> %s" % ("", hits))

# ----------------------------------------------------------------------------
# (3) cross-family identities discovered by scan: equal sum-sequences
# ----------------------------------------------------------------------------
print("\n" + "=" * 100)
print("(3) CROSS-FAMILY COINCIDENCES (identical first-10 sum sequences)")
print("=" * 100)
sigs = {}
for name, g, finite in FAMS:
    for t in ([0] if finite else []) + [1, 2, 3]:
        sums = []
        for m in range(12):
            row = sheared_row(g, t, m)
            if row: sums.append(sum(row))
        key = tuple(sums[:10])
        if len(key) >= 8: sigs.setdefault(key, []).append("%s.t%d" % (name, t))
dups = {k: v for k, v in sigs.items() if len(v) > 1}
if dups:
    for k, v in dups.items():
        print("  MATCH %s : %s..." % (" == ".join(v), list(k[:6])))
else:
    print("  none at first-10 exactness")

# ----------------------------------------------------------------------------
# (4) the q-axis and the relational reading, deepened: [n,2]_q for q = 1..5
#     -- "n in the relational sense" over richer ground structures
# ----------------------------------------------------------------------------
print("\n" + "=" * 100)
print("(4) THE RELATIONAL AXIS DEEPENED: [n choose 2]_q  (q-deformed 'edges of K_n')")
print("=" * 100)
for q1 in range(0, 5):
    col = [fam_qbinom(q1, n) for n in range(2, 11)]
    print("  q=%d: %s  -> %s" % (1 + q1, col, identify(col)))
print("  q->1 is SET relations (triangular); q=2 is F2-LINEAR relations (2-subspaces);")
print("  the q-axis is the fourth continuation: vary the GROUND FIELD of 'relation'.")

print("\nDONE.")
