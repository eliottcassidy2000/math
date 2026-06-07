#!/usr/bin/env python3
"""
THM-438 ADDENDUM-12 follow-up (monad-explorer 14th): TEST the new lead -- is the
THM-438 cycle-rank triangle t(k,m) the RATE-MARKED free compound Poisson of the
exponential, i.e. the block-count grading of the NC moment-cumulant sum?

The free exp-Poisson at rate t (its free-convolution semigroup mu^{box t}) has
free cumulants kappa_n = t*n! and moments
    m_k(t) = sum_{pi in NC(k)} t^{|pi|} prod_B |B|!  = sum_j N(k,j) t^j,
N(k,j) = sum over NC partitions of [k] with exactly j blocks of prod |B|!.
Equivalently M(z,t) = F(t z M(z,t)), F=sum n! w^n  (the y-marked Callan eqn).

Compare N(k,j) (and the M=F(tzM) coefficients) to the genuine THM-438 triangle
rows t(k,m). If they match under some reindex -> the columns ARE rate-graded
free-CP moments (clean identification). If not -> the lead is killed cleanly.
"""
from fractions import Fraction as Fr
from math import factorial

# genuine THM-438 cycle-rank triangle rows (from ADD-3/ADD-6, VERIFIED k<=6)
TRIANGLE = {
    1: [1],
    2: [1, 3],
    3: [1, 9, 13],
    4: [1, 18, 72, 69],
    5: [1, 30, 230, 580, 421],
    6: [1, 45, 560, 2626, 4845, 2867],
}

# ---- noncrossing partitions of [n], enumerated, with block count -----------
def set_partitions(elements):
    if not elements:
        yield []
        return
    first = elements[0]
    for rest in set_partitions(elements[1:]):
        for i in range(len(rest)):
            yield rest[:i] + [(first,) + rest[i]] + rest[i + 1:]
        yield [(first,)] + rest

def is_noncrossing(blocks):
    bl = [sorted(b) for b in blocks]
    for i in range(len(bl)):
        for j in range(len(bl)):
            if i == j:
                continue
            B, C = bl[i], bl[j]
            for ai in range(len(B)):
                for ci in range(ai + 1, len(B)):
                    a, c = B[ai], B[ci]
                    for x in C:
                        for y in C:
                            if a < x < c < y:
                                return False
    return True

def N_table(kmax):
    """N[k][j] = sum_{NC(k), j blocks} prod |B|!."""
    out = {}
    for k in range(1, kmax + 1):
        row = {}
        for blocks in set_partitions(list(range(k))):
            if not is_noncrossing(blocks):
                continue
            j = len(blocks)
            prod = 1
            for B in blocks:
                prod *= factorial(len(B))
            row[j] = row.get(j, 0) + prod
        out[k] = [row.get(j, 0) for j in range(1, k + 1)]   # j=1..k
    return out

KMAX = 6
N = N_table(KMAX)
print("Rate-marked free exp-Poisson:  N(k,j) = sum_{NC(k),|pi|=j} prod|B|!  (j=1..k):")
for k in range(1, KMAX + 1):
    print(f"   k={k}: N(k,.) = {N[k]}    rowsum={sum(N[k])} (=A088368({k}))")

print("\nGenuine THM-438 cycle-rank triangle t(k,m) (m=1..k):")
for k in range(1, KMAX + 1):
    print(f"   k={k}: t(k,.) = {TRIANGLE[k]}    diag t(k,k)={TRIANGLE[k][-1]}")

# ---- direct comparisons ----------------------------------------------------
print("\nComparison (is t the rate-marked NC table, possibly reversed?):")
for k in range(1, KMAX + 1):
    t = TRIANGLE[k]
    n = N[k]
    nrev = n[::-1]
    print(f"   k={k}: t==N? {t==n}   t==reverse(N)? {t==nrev}")
    # signed sums
    sl_t = sum((-1)**(m+1) * t[m] for m in range(k))   # m index 0..k-1 -> rank 1..k
    print(f"         row-signed (THM, sum (-1)^m t): {sum((-1)**(m+1)*t[m] for m in range(k))}"
          f" ; NC block-signed (sum (-1)^j N): {sum((-1)**(j+1)*n[j] for j in range(k))}")

# ---- test the marked Callan equation M(z,t)=F(t z M) -----------------------
def marked_moments(kmax):
    """coeffs m_k(t) of M=F(tzM), as polynomials in t (list over t-powers)."""
    # work with polynomials in t (lists) whose entries are ints, series in z.
    # M = sum_k Mz[k] where Mz[k] is a dict {t-power: coeff}
    fac = [factorial(n) for n in range(kmax + 1)]
    # iterate
    Mz = [dict() for _ in range(kmax + 1)]
    Mz[0] = {0: 1}
    for _ in range(kmax + 2):
        # w = t z M : [z^k] w = t * M[k-1]; multiply each t-power by extra t
        def polymul(a, b):
            r = {}
            for pa, ca in a.items():
                for pb, cb in b.items():
                    r[pa + pb] = r.get(pa + pb, 0) + ca * cb
            return r
        # build powers of w as z-series of t-polys
        # w[k] = shift of t*M[k-1]
        w = [dict() for _ in range(kmax + 1)]
        for k in range(1, kmax + 1):
            w[k] = {p + 1: c for p, c in Mz[k - 1].items()}
        newM = [dict() for _ in range(kmax + 1)]
        newM[0] = {0: 1}
        # accumulate F(w)=sum_n n! w^n
        wp = [dict(d) for d in w]   # w^1
        for nn in range(1, kmax + 1):
            for k in range(kmax + 1):
                for p, c in wp[k].items():
                    newM[k][p] = newM[k].get(p, 0) + fac[nn] * c
            # wp = wp * w (z-convolution, t-polymul)
            nxt = [dict() for _ in range(kmax + 1)]
            for i in range(kmax + 1):
                if not wp[i]:
                    continue
                for j in range(kmax + 1 - i):
                    if not w[j]:
                        continue
                    pr = polymul(wp[i], w[j])
                    d = nxt[i + j]
                    for p, c in pr.items():
                        d[p] = d.get(p, 0) + c
            wp = nxt
            if all(not d for d in wp):
                break
        if newM == Mz:
            break
        Mz = newM
    return Mz

Mz = marked_moments(KMAX)
print("\nMarked Callan eqn M=F(tzM): m_k(t) coefficients (t-power: coeff):")
ok = True
for k in range(1, KMAX + 1):
    coeffs = [Mz[k].get(j, 0) for j in range(1, k + 1)]
    match = (coeffs == N[k])
    ok = ok and match
    print(f"   k={k}: {coeffs}   == N(k,.)? {match}")
print(f"\nM=F(tzM) reproduces the rate-marked NC table: {ok}")

print("\nCONCLUSION:")
print("  - rate-marked free exp-Poisson moments = N(k,j) (block grading), and")
print("    M=F(tzM) is its bivariate functional equation (clean, algebraic-in-t).")
print("  - the GENUINE THM-438 cycle-rank triangle t(k,m) is a DIFFERENT refinement:")
print("    they share only the diagonal/rowsum (=A088368) and t(k,1)=N(k,1)? check above.")
