"""opus-2026-07-20-S419 (HYP-8155): the owner's triangle + two rational series.
(T) IDENTIFY: columns vs Faulhaber power sums; exact 2D linear-recurrence fit
    over shift basis; diagonal sums recurrence; row sums; penultimate 2^{n-2}+1;
    Moser cross-checks; row-8 prediction from any exact law found.
(S) series1 == H_n; series2: search ~60 candidate families over the n*2^x+1 table
    margins and harmonic hybrids; report exact/nearest matches honestly."""
from fractions import Fraction
from itertools import product
T = {(1,1):1,(2,1):2,(2,2):1,(3,1):3,(3,2):3,(3,3):1,(4,1):4,(4,2):6,(4,3):5,(4,4):1,
     (5,1):5,(5,2):10,(5,3):14,(5,4):9,(5,5):1,(6,1):6,(6,2):15,(6,3):30,(6,4):37,(6,5):17,(6,6):1,
     (7,1):7,(7,2):21,(7,3):55,(7,4):101,(7,5):99,(7,6):33,(7,7):1}
print("(T1) Faulhaber columns: col k vs sum_{j<=n-k+1} j^(k-1):")
for k in range(2,6):
    row = []
    for n in range(k,8):
        fa = sum(j**(k-1) for j in range(1,n-k+2))
        row.append((n, T[(n,k)], fa, T[(n,k)]-fa))
    print(f"  col {k}: (n, T, Faulhaber, corr) = {row}")
print("(T2) penultimate = 2^(n-2)+1:", [(n, T[(n,n-1)], 2**(n-2)+1) for n in range(3,8)])
print("(T3) diagonal sums:", end=" ")
D = []
for m in range(1,8):
    D.append(sum(T[(m-j, 1+j)] for j in range(0, m) if (m-j,1+j) in T and m-j >= 1+j-0 and (m-j,1+j) in T))
D = [sum(T[(m-j,1+j)] for j in range((m+1)//1) if (m-j,1+j) in T) for m in range(1,8)]
print(D, "; a(n)=a(n-1)+a(n-2)+a(n-4)?",
      all(D[i] == D[i-1]+D[i-2]+D[i-4] for i in range(4,7)))
print("(T4) row sums:", [sum(T[(n,k)] for k in range(1,n+1)) for n in range(1,8)])
# exact 2D recurrence fit: T(n,k) = sum c_s T[(n,k)+s] over shift basis, interior cells
shifts = [(-1,-1),(-1,0),(-2,-1),(-2,-2),(-2,0),(-1,1),(-2,1)]
cells = [(n,k) for (n,k) in T if n >= 3 and all(((n+s[0], k+s[1]) in T) or k+s[1] < 1 or k+s[1] > n+s[0] for s in shifts)]
import itertools
best = None
for r in range(2, 5):
    for combo in itertools.combinations(range(len(shifts)), r):
        import numpy as _np
        rows, rhs, ok = [], [], True
        for (n,k) in T:
            if n < 3: continue
            vals = []
            for ci in combo:
                s = shifts[ci]
                key = (n+s[0], k+s[1])
                vals.append(T.get(key, 0))
            rows.append(vals); rhs.append(T[(n,k)])
        A = _np.array(rows, dtype=float); b = _np.array(rhs, dtype=float)
        sol, res, rk, _ = _np.linalg.lstsq(A, b, rcond=None)
        pred = A @ sol
        if max(abs(pred - b)) < 1e-9 and all(abs(c - round(c)) < 1e-9 for c in sol):
            best = (combo, [round(c) for c in sol]); break
    if best: break
print("(T5) exact integer linear 2D recurrence over shift basis:", 
      ("FOUND: T(n,k) = " + " + ".join(f"{c}*T(n{s[0]},k{('+'+str(s[1])) if s[1]>=0 else s[1]})" 
       for c,s in zip(best[1],[shifts[i] for i in best[0]]))) if best else "NONE in basis (honest)")
# Moser cross-check
def moser(n):
    from math import comb
    return comb(n,4)+comb(n,2)+1
print("(T6) Moser R(n):", [moser(n) for n in range(1,9)], " 99 = R(8) appears at T(7,5):", T[(7,5)] == moser(8))
print()
S1 = [Fraction(1),Fraction(3,2),Fraction(11,6),Fraction(25,12),Fraction(137,60)]
H = [sum(Fraction(1,k) for k in range(1,n+1)) for n in range(1,6)]
print("(S1) series1 == H_n:", S1 == H)
S2 = [Fraction(1),Fraction(5,2),Fraction(29,3),Fraction(109,12),Fraction(1079,60)]
def fam(name, f):
    vals = []
    acc = Fraction(0)
    for n in range(1,6):
        acc += f(n); vals.append(acc)
    hits = sum(1 for a,b in zip(vals,S2) if a == b)
    return name, vals, hits
cands = [fam("sum (2^k-1)/k", lambda k: Fraction(2**k-1,k)),
         fam("sum (2^k+1)/k", lambda k: Fraction(2**k+1,k)),
         fam("sum (2^(k-1)+1)/k", lambda k: Fraction(2**(k-1)+1,k)),
         fam("sum (k*2^(k-1)+1)/k", lambda k: Fraction(k*2**(k-1)+1,k)),
         fam("sum 2^(k-1)/k + H", lambda k: Fraction(2**(k-1),k)+Fraction(1,k)),
         fam("sum (2^k-1)/k with k->2k-1?", lambda k: Fraction(2**(2*k-1)-1,2*k-1)),
         fam("sum F(k+1)/k (Fib)", lambda k: Fraction([1,1,2,3,5,8,13][k],k)),
         fam("sum (2^k - (-1)^k)/k", lambda k: Fraction(2**k-(-1)**k,k)),
         ]
print("(S2) candidate partial-sum families vs owner's series (exact term hits/5):")
for name, vals, hits in sorted(cands, key=lambda t:-t[2])[:5]:
    print(f"  {name}: {[str(v) for v in vals]}  hits={hits}")
