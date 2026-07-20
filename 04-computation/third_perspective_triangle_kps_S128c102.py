#!/usr/bin/env python3
"""
kind-pasteur-2026-07-20-S128c102 -- HYP-8165: the owner's THIRD-PERSPECTIVE
TRIANGLE.  OEIS: triangle flat, row sums, diag sums, series-2 numerators ALL
ABSENT (curl verified) -- identification from within.

 (1) exact deviation table vs power sums Sum_{j=1}^{m} j^(k-1), m = n-k+1
 (2) exhaustive small-support linear recurrence fitter (constant + k-linear
     + n-linear coefficients)
 (3) row sums, alternating row sums, diagonal (NE) sums + Fibonacci probes
 (4) series-2 mega-fitter with leave-one-out typo detection
 (5) the Proth/repo crosswalk: 2n+1 (H-spectrum odds), 2^x+1 (hypotenuse law,
     Fermat rungs 3,5,17,257), Moser positions, penultimate column
"""
import sympy as sp
from sympy import Rational, symbols, binomial, fibonacci
from itertools import product as iproduct
from fractions import Fraction as Fr

T = {
 (1,1): 1,
 (2,1): 2, (2,2): 1,
 (3,1): 3, (3,2): 3, (3,3): 1,
 (4,1): 4, (4,2): 6, (4,3): 5, (4,4): 1,
 (5,1): 5, (5,2): 10, (5,3): 14, (5,4): 9, (5,5): 1,
 (6,1): 6, (6,2): 15, (6,3): 30, (6,4): 37, (6,5): 17, (6,6): 1,
 (7,1): 7, (7,2): 21, (7,3): 55, (7,4): 101, (7,5): 99, (7,6): 33, (7,7): 1,
}
NMAX = 7

print("== (1) deviations vs power sums ==", flush=True)
for (n,k), v in sorted(T.items()):
    m = n - k + 1
    ps = sum(j**(k-1) for j in range(1, m+1))
    d = v - ps
    if d != 0:
        print(f"  T({n},{k}) = {v} = Sum_(j<={m}) j^{k-1} {'+' if d>0 else '-'} {abs(d)}   [m={m}]", flush=True)
print("  (all other entries: deviation 0 -- cols 1,2,3 pure; penultimate = 1+2^(n-2) = the 2-term power sum)", flush=True)

print("\n== (2) recurrence fitter ==", flush=True)
# candidate supports: shifts (dn, dk) with dn in 1..2, and coefficient forms 1, k, n, n-k
shifts = [(1,0),(1,1),(2,0),(2,1),(2,2)]
coefkinds = ['1','k','n','m']  # m = n-k+1
import itertools
def coefval(kind, n, k):
    return {'1':1,'k':k,'n':n,'m':n-k+1}[kind]
found = []
# try supports of size 2 and 3 with unknown rational constant coefficients times a coefficient-kind
for r in (2, 3):
    for supp in itertools.combinations(shifts, r):
        for kinds in iproduct(coefkinds, repeat=r):
            # solve linear system for constants c_i: T(n,k) = sum c_i * kind_i(n,k) * T(n-dn_i, k-dk_i)
            rows = []
            rhs = []
            for (n,k), v in T.items():
                ok = True
                row = []
                for (dn,dk), kd in zip(supp, kinds):
                    nn, kk = n-dn, k-dk
                    if kk < 1 or kk > nn or nn < 1:
                        val = 0
                    else:
                        val = T[(nn,kk)]
                    row.append(Fr(coefval(kd, n, k)) * val)
                if n <= 2: continue   # skip boundary rows
                rows.append(row); rhs.append(Fr(v))
            # least-squares-free exact solve: use first r independent rows then verify all
            import sympy as sp2
            A = sp2.Matrix([[sp2.Rational(x.numerator, x.denominator) for x in row] for row in rows])
            bv = sp2.Matrix([sp2.Rational(x.numerator, x.denominator) for x in rhs])
            try:
                sol, params = A.gauss_jordan_solve(A[:,:], bv)
            except Exception:
                continue
            if params.shape[0] > 0:  # underdetermined: skip
                continue
            resid = A*sol - bv
            if all(e == 0 for e in resid):
                found.append((supp, kinds, list(sol)))
                print(f"  EXACT RECURRENCE: T(n,k) = " + " + ".join(
                    f"({c})*{kd}*T(n-{dn},k-{dk})" for (dn,dk), kd, c in zip(supp, kinds, sol)), flush=True)
if not found:
    print("  no constant/k/n/m-coefficient linear recurrence on these supports (sizes 2-3)", flush=True)

print("\n== (3) sums and Fibonacci probes ==", flush=True)
rs = [sum(T[(n,k)] for k in range(1, n+1)) for n in range(1, NMAX+1)]
als = [sum((-1)**(k+1)*T[(n,k)] for k in range(1, n+1)) for n in range(1, NMAX+1)]
print(f"  row sums: {rs}", flush=True)
print(f"  alternating row sums: {als}", flush=True)
diags = []
for m in range(1, NMAX+1):
    sdi = 0
    j = 0
    while True:
        n = m - j; k = j + 1
        if n < 1 or k > n: break
        sdi += T[(n,k)]
        j += 1
    diags.append(sdi)
print(f"  NE-diagonal sums (Pascal->Fibonacci analog): {diags}", flush=True)
# Fibonacci residual: diag - fib?
fibs = [sp.fibonacci(i) for i in range(1, NMAX+1)]
print(f"  Fibonacci:                                    {fibs}", flush=True)
print(f"  diag - Fibonacci: {[d-f for d,f in zip(diags, fibs)]}", flush=True)
# weighted diagonal (Moser-style even strides)
print(f"  2^n:      {[2**(n-1) for n in range(1, NMAX+1)]}", flush=True)
print(f"  Moser:    {[1 + sp.binomial(n-1,2) + sp.binomial(n-1,4) for n in range(1, NMAX+1)]}", flush=True)
print(f"  rowsum - 2^(n-1): {[r - 2**(n-1) for r, n in zip(rs, range(1, NMAX+1))]}", flush=True)
print(f"  rowsum - Moser:   {[r - (1 + sp.binomial(n-1,2) + sp.binomial(n-1,4)) for r, n in zip(rs, range(1, NMAX+1))]}", flush=True)

print("\n== (4) series-2 mega-fitter (owner: 1, 5/2, 29/3, 109/12, 1079/60) ==", flush=True)
targ = [Fr(1), Fr(5,2), Fr(29,3), Fr(109,12), Fr(1079,60)]
targ_alt = [Fr(1), Fr(5,2), Fr(29,6), Fr(109,12), Fr(1079,60)]  # 29/6 typo theory
def check_terms(name, termfn, tgt):
    acc = Fr(0); good = 0
    detail = []
    for k in range(1, 6):
        acc += termfn(k)
        detail.append(acc)
        if acc == tgt[k-1]: good += 1
    return good, detail
cands = {
 "Sum (2^k-1)/k [= Sum C(n,k)/k]": lambda k: Fr(2**k - 1, k),
 "Sum (2^(k-1)+1)/k": lambda k: Fr(2**(k-1) + 1, k),
 "Sum 2^(k-1)/k + H": lambda k: Fr(2**(k-1), k) + Fr(1, k),
 "Sum (2^k+(-1)^k)/k": lambda k: Fr(2**k + (-1)**k, k),
 "Sum (k*2^(k-1)+1)/k^2?": lambda k: Fr(k*2**(k-1) + 1, k*k),
 "Sum penult(k+1)/k = (2^(k-1)+1)/k": lambda k: Fr(2**(k-1)+1, k),
 "Sum (2^k - k)/k": lambda k: Fr(2**k - k, k),
 "Sum C(2k,k)/(k+1)/k?": lambda k: Fr(int(sp.binomial(2*k,k)), (k+1)*k),
}
for tgt, lbl in ((targ, "as-given"), (targ_alt, "29/6-typo")):
    print(f"  target ({lbl}): {[str(t) for t in tgt]}", flush=True)
    best = []
    for name, f in cands.items():
        good, detail = check_terms(name, f, tgt)
        if good >= 3:
            print(f"    {name}: matches {good}/5; partials {[str(d) for d in detail]}", flush=True)
# direct term extraction
for tgt, lbl in ((targ, "as-given"), (targ_alt, "29/6-typo")):
    terms = [tgt[0]] + [tgt[i]-tgt[i-1] for i in range(1, 5)]
    print(f"  {lbl}: term differences: {[str(t) for t in terms]}", flush=True)

print("\n== (5) the Proth/repo crosswalk ==", flush=True)
print("  penultimate column: " + str([T[(n,n-1)] for n in range(2, 8)]) + " = 1+2^(n-2) = THE HYPOTENUSE LAW / transitive-SC-neighbor H", flush=True)
print("  prime entries: 2, 3, 5, 17 -- exactly the THM-871 Fermat rungs (3, 5, 17, 257 continue); 9 = 3^2, 33 = 3*11 composite", flush=True)
print("  H-spectrum (repo): Hamiltonian-path counts = odd numbers 2n+1 minus {7, 21} -- the owner's (1,n) row", flush=True)
print("  n*2^x+1 = Proth numbers: the table interpolates the H-spectrum frame (x=1) and the hypotenuse tower (n=1)", flush=True)
print("  T(7,2) = 21 and T(3,1) = 3: the two FORBIDDEN H-values {7, 21}: 7 = row-3 sum = 2^3-1, 21 = C(7,2) = T(7,2)", flush=True)
print("\nDONE.", flush=True)
