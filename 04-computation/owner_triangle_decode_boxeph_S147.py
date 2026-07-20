#!/usr/bin/env python3
"""
owner_triangle_decode_boxeph_S147.py  (HYP-8165)

DECODE the owner triangle
  {1},{2,1},{3,3,1},{4,6,5,1},{5,10,14,9,1},{6,15,30,37,17,1},{7,21,55,101,99,33,1}
and verify/mine every claimed connection: triangular/pyramidal columns, 2^x+1
diagonal, Fibonacci, powers of 2, Moser, the harmonic series H_n, the second
rational series, and the (2n+1, 2^x+1, n*2^x+1) grid.  Exact arithmetic.

boxeph-2026-07-20-S147.
"""

from fractions import Fraction as Fr
from itertools import product

T = {1: [1], 2: [2,1], 3: [3,3,1], 4: [4,6,5,1], 5: [5,10,14,9,1],
     6: [6,15,30,37,17,1], 7: [7,21,55,101,99,33,1]}
def t(m, j):
    if m in T and 1 <= j <= m: return T[m][j-1]
    return 0

print("=" * 92)
print("(A) RECURRENCE FIT: stencil c1*T(m-1,j-1) + c2*T(m-1,j) + c3*T(m-2,j-1) + ...")
print("=" * 92)
# model: T(m,j) = sum over stencil of (a_s + b_s*j + c_s*(m-j)) * T(m-ds, j-es)
stencils = [
    [(1,1),(1,0)],
    [(1,1),(1,0),(2,1)],
    [(1,1),(1,0),(2,2)],
    [(1,1),(1,0),(2,1),(2,2)],
    [(1,1),(1,0),(2,0),(2,1),(2,2)],
]
best = None
for stn in stencils:
    for use_j in (False, True):
        for use_k in (False, True):   # k = m - j
            nun = len(stn) * (1 + use_j + use_k)
            rows, rhs = [], []
            for m in range(2, 8):
                for j in range(1, m+1):
                    row = []
                    for (ds, es) in stn:
                        v = t(m-ds, j-es)
                        row.append(Fr(v))
                        if use_j: row.append(Fr(v*j))
                        if use_k: row.append(Fr(v*(m-j)))
                    rows.append(row); rhs.append(Fr(t(m,j)))
            # exact least-squares-free solve: Gaussian on the system; check consistency
            A = [r[:] + [rhs[i]] for i, r in enumerate(rows)]
            nc = nun; rk = 0; piv = []
            for c in range(nc):
                p = None
                for r in range(rk, len(A)):
                    if A[r][c] != 0: p = r; break
                if p is None: continue
                A[rk], A[p] = A[p], A[rk]
                pr = A[rk]; inv = 1/pr[c]
                for r in range(len(A)):
                    if r != rk and A[r][c] != 0:
                        f = A[r][c]*inv
                        for cc in range(nc+1): A[r][cc] -= f*pr[cc]
                piv.append(c); rk += 1
            ok = all(all(A[r][c] == 0 for c in range(nc)) is False or A[r][nc] == 0 for r in range(len(A)))
            ok = True
            for r in range(len(A)):
                if all(A[r][c] == 0 for c in range(nc)) and A[r][nc] != 0: ok = False
            if ok:
                sol = [Fr(0)]*nc
                for i, c in enumerate(piv): sol[c] = A[i][nc]/A[i][c]
                # verify perfectly
                good = True
                for m in range(2, 8):
                    for j in range(1, m+1):
                        acc = Fr(0); ci = 0
                        for (ds, es) in stn:
                            v = t(m-ds, j-es)
                            acc += sol[ci]*v; ci += 1
                            if use_j: acc += sol[ci]*v*j; ci += 1
                            if use_k: acc += sol[ci]*v*(m-j); ci += 1
                        if acc != t(m, j): good = False
                if good and (best is None or nun < best[0]):
                    best = (nun, stn, use_j, use_k, sol)
if best:
    nun, stn, use_j, use_k, sol = best
    print("EXACT RULE FOUND: stencil %s, coeffs (const%s%s) = %s" %
          (stn, "/xj" if use_j else "", "/xk" if use_k else "", [str(s) for s in sol]))
    def extend(rule, upto=12):
        TT = {m: T[m][:] for m in T}
        for m in range(8, upto+1):
            row = []
            for j in range(1, m+1):
                acc = Fr(0); ci = 0
                for (ds, es) in stn:
                    v = TT.get(m-ds, [0]*(m))[j-es-1] if 1 <= j-es <= m-ds else 0
                    acc += sol[ci]*v; ci += 1
                    if use_j: acc += sol[ci]*v*j; ci += 1
                    if use_k: acc += sol[ci]*v*(m-j); ci += 1
                row.append(int(acc))
            TT[m] = row
        return TT
    TT = extend(best)
    print("row 8:", TT[8])
    print("row 9:", TT[9])
else:
    print("no exact rule among local stencils; trying SHIFTED-BINOMIAL basis...")
    from math import comb
    basis = [(a, b) for a in range(0, 4) for b in range(0, 4)]
    A = []
    for m in range(1, 8):
        for j in range(1, m+1):
            A.append([Fr(comb(max(m-a,0), j-b)) if 0 <= j-b <= m-a else Fr(0)
                      for (a, b) in basis] + [Fr(t(m,j))])
    nc = len(basis); rk = 0; piv = []
    for c in range(nc):
        pv = None
        for r in range(rk, len(A)):
            if A[r][c] != 0: pv = r; break
        if pv is None: continue
        A[rk], A[pv] = A[pv], A[rk]
        pr = A[rk]; inv = 1/pr[c]
        for r in range(len(A)):
            if r != rk and A[r][c] != 0:
                f = A[r][c]*inv
                for cc in range(nc+1): A[r][cc] -= f*pr[cc]
        piv.append(c); rk += 1
    bad = any(all(A[r][c] == 0 for c in range(nc)) and A[r][nc] != 0 for r in range(len(A)))
    if not bad:
        sol = [Fr(0)]*nc
        for i, c in enumerate(piv): sol[c] = A[i][nc]/A[i][c]
        nz = [(basis[i], sol[i]) for i in range(nc) if sol[i]]
        print("EXACT BINOMIAL FORM: T(m,j) = " +
              " + ".join("%s*C(m-%d,j-%d)" % (c, a, b) for (a,b), c in nz))
        from math import comb as CB
        def tval(m, j):
            return sum(int(c)*CB(m-a, j-b) if 0 <= j-b <= m-a else 0 for (a,b), c in nz)
        ok = all(tval(m,j) == t(m,j) for m in range(1,8) for j in range(1,m+1))
        print("verified on all 28 entries:", ok)
        if ok:
            for m in (8, 9):
                T[m] = [tval(m, j) for j in range(1, m+1)]
                print("row %d:" % m, T[m])
    else:
        print("binomial basis (16 shifts) also insufficient -- REPORT AS OPEN: the")
        print("triangle is not a finite shifted-binomial combination at these shifts.")
    TT = {m: T[m][:] for m in T}

print("\n" + "=" * 92)
print("(B) STRUCTURE CHECKS (exact)")
print("=" * 92)
tri = [m*(m+1)//2 for m in range(0, 10)]
pyr = [m*(m+1)*(2*m+1)//6 for m in range(0, 10)]
print("col1 = n:", [t(m,1) for m in range(1,8)])
print("col2 = triangular:", [t(m,2) for m in range(2,8)], "==", tri[1:7])
print("col3 = square pyramidal:", [t(m,3) for m in range(3,8)], "==", pyr[1:6])
print("diag (m,m-1) = 2^(m-2)+1:", [t(m,m-1) for m in range(2,8)],
      "==", [2**(m-2)+1 for m in range(2,8)])
rows_sums = {m: sum(T[m]) for m in T}
print("row sums:", [rows_sums[m] for m in range(1,8)])
alt = {m: sum((-1)**(j-1)*t(m,j) for j in range(1,m+1)) for m in T}
print("alternating row sums:", [alt[m] for m in range(1,8)])
# shallow diagonal sums (Fibonacci-analogue reading)
sh = []
for s in range(1, 10):
    v = sum(t(s-i, 1+i) for i in range(0, s))
    sh.append(v)
print("shallow-diagonal sums:", sh, " (Pascal's would give Fibonacci)")
fib = [1,1,2,3,5,8,13,21,34,55]
print("Fibonacci entries appearing in the triangle:",
      sorted({t(m,j) for m in T for j in range(1,m+1)} & set(fib)))
# Moser reading: partial row sums truncated at column 3 (1 + C(n,2) + C(n,4)-analogue)
print("partial row sums cols 1..3:", [sum(T[m][:3]) for m in range(1,8)],
      " vs Moser 1,2,4,8,16,31,57")

print("\n" + "=" * 92)
print("(C) THE RATIONAL SERIES")
print("=" * 92)
H = [sum(Fr(1,k) for k in range(1,n+1)) for n in range(1,8)]
print("series1 owner: 1, 3/2, 11/6, 25/12, 137/60  == H_n:", [str(h) for h in H[:5]])
owner2 = [Fr(1), Fr(5,2), Fr(29,3), Fr(109,12), Fr(1079,60)]
cands = {
 "sum rowsum(k)/k":      [sum(Fr(rows_sums[k],k) for k in range(1,n+1)) for n in range(1,6)],
 "sum T(n,j)/j (row n)": [sum(Fr(t(n,j),j) for j in range(1,n+1)) for n in range(1,6)],
 "sum T(n,j)/(n+1-j)":   [sum(Fr(t(n,j),n+1-j) for j in range(1,n+1)) for n in range(1,6)],
 "sum rowsum(k)/2^(k-1)":[sum(Fr(rows_sums[k],2**(k-1)) for k in range(1,n+1)) for n in range(1,6)],
 "H_n * n-ish 2H_n-1":   [2*H[n-1]-Fr(1) for n in range(1,6)],
 "sum (2^k-1)/k":        [sum(Fr(2**k-1,k) for k in range(1,n+1)) for n in range(1,6)],
 "sum t(k,k-1)+.../k":   [sum(Fr(t(k+1,k),k) for k in range(1,n+1)) for n in range(1,6)],
}
for name, vals in cands.items():
    hits = sum(1 for a, b in zip(vals, owner2) if a == b)
    print("  %-24s -> %s  matches owner2: %d/5" % (name, [str(v) for v in vals], hits))
print("  owner2 =", [str(v) for v in owner2], " (NB 29/3 vs candidates' 29/6: possible typo)")

print("\n" + "=" * 92)
print("(D) THE (x,n) GRID n*2^x+1 -- Proth/Sierpinski face + owner's boundary spec")
print("=" * 92)
print("     n*2^x+1 grid (rows n=0..5, cols x=0..5):")
for n in range(0, 6):
    print("  n=%d:" % n, [n*2**x+1 for x in range(0, 6)])
print("  owner boundary spec: (1,n)=2n+1 OK (x=1 col); (x,1)=2^x+1 OK (n=1 row);")
print("  (0,n)=n+1 (owner says n) and (x,0)=1 OK -- the owner's chart is the n*2^x+1")
print("  grid up to the (0,n) off-by-one, i.e. the PROTH/SIERPINSKI table: 2n+1 =")
print("  our pigeonhole moduli; 2^x+1 = Fermat/gate primes; the x-direction = the")
print("  x2-tower of the S-T gates; Sierpinski numbers = rows with NO primes.")
print("DONE.")
