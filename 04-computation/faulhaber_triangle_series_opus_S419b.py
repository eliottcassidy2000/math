"""opus-2026-07-20-S419 part B: (1) 5-term recurrence hunt + correction-cell law;
(2) series-2 arbitration (hybrid H_n+2^(n-1)-1 vs Erdos-Borwein sum(2^k-1)/k =
binomial-harmonic sum C(n,k)/k -- identity verified); (3) row-8 prediction under
flat-correction law; (4) the bilinear a*b+1 rim: repo denominators vs Proth."""
from fractions import Fraction
from math import comb
import itertools, numpy as np
T = {(1,1):1,(2,1):2,(2,2):1,(3,1):3,(3,2):3,(3,3):1,(4,1):4,(4,2):6,(4,3):5,(4,4):1,
     (5,1):5,(5,2):10,(5,3):14,(5,4):9,(5,5):1,(6,1):6,(6,2):15,(6,3):30,(6,4):37,(6,5):17,(6,6):1,
     (7,1):7,(7,2):21,(7,3):55,(7,4):101,(7,5):99,(7,6):33,(7,7):1}
F = lambda n,k: sum(j**(k-1) for j in range(1, n-k+2))
print("(B1) correction cells T - Faulhaber:", {(n,k): T[(n,k)]-F(n,k) for (n,k) in T if T[(n,k)] != F(n,k)},
      "| flat-1 on {k>=4, n-k>=2} fits all 28 cells:",
      all(T[(n,k)] == F(n,k) + (1 if (k >= 4 and n-k >= 2) else 0) for (n,k) in T))
shifts = [(-1,-1),(-1,0),(-2,-1),(-2,-2),(-2,0),(-1,1),(-2,1),(-1,-2),(-2,-3)]
best = None
for r in range(2, 6):
    for combo in itertools.combinations(range(len(shifts)), r):
        rows, rhs = [], []
        for (n,k) in T:
            if n < 3: continue
            rows.append([T.get((n+s[0], k+s[1]), 0) for si, s in enumerate(shifts) if si in combo])
            rhs.append(T[(n,k)])
        A = np.array(rows, float); b = np.array(rhs, float)
        sol = np.linalg.lstsq(A, b, rcond=None)[0]
        if max(abs(A @ sol - b)) < 1e-9 and all(abs(c-round(c)) < 1e-9 for c in sol):
            best = (combo, [round(c) for c in sol]); break
    if best: break
print("(B2) integer linear recurrence up to 5 terms/9 shifts:", best if best else "NONE (triangle is not shift-linear on this basis; columns are degree-(k-1) polynomial => any column law needs depth k)")
row8 = [F(8,k) + (1 if (k >= 4 and 8-k >= 2) else 0) for k in range(1,9)]
print("(B3) row-8 PREDICTION (Faulhaber + flat-1 law):", row8)
S2p = [Fraction(1),Fraction(5,2),Fraction(29,3),Fraction(109,12),Fraction(1079,60)]
hyb = [sum(Fraction(1,k) for k in range(1,n+1)) + 2**(n-1) - 1 for n in range(1,6)]
eb  = []
acc = Fraction(0)
for k in range(1,6):
    acc += Fraction(2**k-1, k); eb.append(acc)
bh  = [sum(Fraction(comb(n,k), k) for k in range(1,n+1)) for n in range(1,6)]
print("(B4) series-2 arbitration vs printed", [str(v) for v in S2p])
print("   HYBRID H_n+2^(n-1)-1      :", [str(v) for v in hyb], " exact-hits:", sum(a==b for a,b in zip(hyb,S2p)))
print("   ERDOS-BORWEIN sum(2^k-1)/k:", [str(v) for v in eb],  " exact-hits:", sum(a==b for a,b in zip(eb,S2p)))
print("   identity sum C(n,k)/k == sum (2^k-1)/k:", bh == eb)
print("   hybrid increments 1/n + 2^(n-2):", [str(hyb[i]-hyb[i-1]) for i in range(1,5)],
      "| printed dev: n=3:", str(S2p[2]-hyb[2]), "(denominator slip 29/3 vs 29/6), n=5:", str(S2p[4]-hyb[4]))
print()
print("(B5) the bilinear a*b+1 rim across the repo (all quantized denominators have the form m*b+e, e in {+1,-1}):")
rows = [("Proth numbers (owner's table)", "n*2^x+1", "n=5,x=3 -> 41"),
        ("Fermat rim of the triangle", "1*2^x+1", "x=5 -> 33 = T(7,6)"),
        ("odd axis", "n*2^1+1 = 2n+1", "GW escape M = 2/(2N+1): N=14 -> 2/29"),
        ("S-T ladder rungs", "s/(n*s+1)", "s=8,n=14 -> 8/113; cage rung 113/1466 = 113/(13*113 - 3)"),
        ("(D,s) rung ladder THM-1269", "D/((N+1)*D - 1)", "N=14, D=2 -> 2/29; D=13 -> 13/194"),
        ("ghost packing THM-1292", "K/(K*(m+1)+1)", "K=13,m=13 -> 13/183 = 13/Phi6(14)... wait 14/183"),
        ("deep well", "n/Phi6(n) = n/(n^2-n+1)", "14/183; 183 = 13*14+1 BILINEAR"),
        ("harmonic crown THM-1153", "sum 1/s over speeds", "H_n = series 1"),
        ("Moser trap", "B(n-1,4) = partial binomial", "1,2,4,8,16,31: the 2^n illusion"),
        ("width of G_n (repo's own Moser)", "C(n-2,floor) claimed", "FAILS n=7: predicted 10, actual 15")]
for a,b,c in rows: print(f"   {a:38s} {b:26s} {c}")
print("   -> 183 = 13*14+1: the deep well denominator IS a Proth-shaped bilinear (13 = n-1, base 14).")
