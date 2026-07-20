#!/usr/bin/env python3
"""
death-star-2026-07-20-S59t (HYP-8175) -- the n*2^x+1 staggerings (sum & product,
slopes 0/1/2), the RELATIONAL reading (triangular = edges of K_n = tournament
arcs = staircase tiles; 2^triangular = graphs/tilings), and the growth-constant
comparison across the whole zoo.
"""
from math import comb, log, isqrt
from fractions import Fraction as Fr
import itertools

def fit(seq,maxord=4):
    seq=list(seq)
    for order in range(1,maxord+1):
        if len(seq)<2*order+1: continue
        for coeffs in itertools.product(range(-3,4),repeat=order):
            if all(c==0 for c in coeffs): continue
            if all(seq[i]==sum(coeffs[j]*seq[i-1-j] for j in range(order)) for i in range(order,len(seq))):
                return coeffs
    return None
def growth(seq):
    seq=[abs(x) for x in seq if x]
    return seq[-1]/seq[-2] if len(seq)>=2 and seq[-2] else None
def charroot(coeffs):
    # dominant real root of x^k - c0 x^{k-1} - ... = 0
    if not coeffs: return None
    x=2.0
    for _ in range(200):
        k=len(coeffs)
        f=x**k - sum(coeffs[j]*x**(k-1-j) for j in range(k))
        df=k*x**(k-1) - sum(coeffs[j]*(k-1-j)*x**(k-2-j) for j in range(k) if k-1-j>=1)
        if df==0: break
        x=x-f/df
    return x

print("=== THE n*2^x+1 GRID staggered (f(x,n)=n*2^x+1) ===")
def f(x,n): return n*2**x+1
# slope-0 (straight grid): row-x sums over n, col-n sums over x
print("  straight grid f(x,n), x=rows 0..5, n=cols 0..6:")
for x in range(6):
    print("   x=%d:"%x,[f(x,n) for n in range(7)])
# antidiagonal slope-s: sum/product over {x+? } ; slope s means row x = m - s*n
def antidiag(m,s,op):
    vals=[f(m-s*n, n) for n in range(0,m+1) if m-s*n>=0]
    if op=="sum": return sum(vals)
    p=1
    for v in vals: p*=v
    return p
for s in (1,2):
    ssum=[antidiag(m,s,"sum") for m in range(0,10)]
    print(f"  slope {s} SUM : {ssum}   rec={fit(ssum)} growth~{growth(ssum):.4f}" if growth(ssum) else "")
for s in (1,2):
    sprod=[antidiag(m,s,"prod") for m in range(0,7)]
    print(f"  slope {s} PROD: {sprod[:6]}")

print("\n=== THE RELATIONAL READING: triangular = the relation itself ===")
T=lambda n: n*(n+1)//2
print("  T_n (triangular):", [T(n) for n in range(9)])
print("  T_{n-1} = C(n,2) = |E(K_n)| = # arcs of a tournament on n vertices:")
print("   ", [(n, comb(n,2)) for n in range(2,9)])
print("  staircase delta_{n-2} tile count = C(n-1,2) = T_{n-2}  (the tiling triangle IS triangular-sized):")
print("   ", [(n, comb(n-1,2)) for n in range(3,9)])
print("  2^{C(n,2)} = 2^triangular = # LABELED GRAPHS on n vertices:")
print("   ", [2**comb(n,2) for n in range(1,7)], "(= 1,2,8,64,1024,32768)")
print("  2^{C(n-1,2)} = # TILINGS (fixed base path) = the tournament tiling hypercube Q_{T_{n-2}}:")
print("   ", [2**comb(n-1,2) for n in range(2,8)])
print("  => x = triangular number is the NATURAL graph/tournament value of the n*2^x+1 exponent:")
print("     n=1, x=T_m gives 2^{T_m}+1 = (#tilings)+1 = the observer-augmented tiling count.")
for m in range(1,6):
    print(f"     2^T_{m}+1 = 2^{T(m)}+1 = {2**T(m)+1}")

print("\n=== GROWTH-CONSTANT LADDER (staggering figurate triangles) ===")
def simplicial(m,c): return comb(m-1+c,c+1) if m>=1 else 0
def pascal(m,c): return comb(m,c)
ladders={
 "Pascal slope0 (rows)":([2**n for n in range(1,10)],"2  (powers of 2)"),
 "Pascal slope1":([sum(comb(m-c,c) for c in range(m+1)) for m in range(1,12)],"phi=1.6180 (Fibonacci)"),
 "Pascal slope2":([sum(comb(m-2*c,c) for c in range(m//2+2)) for m in range(1,12)],"plastic=1.3247 (Narayana/Padovan)"),
 "simplicial slope1":([sum(simplicial(m-c,c) for c in range(m+1)) for m in range(1,10)],"2 (=2^n-1 Mersenne)"),
 "cumulative-Fibonacci":([1,2,4,7,12,20,33,54,88],"phi (Fib partial sums, = simpl.slope2 = k-ary.slope1)"),
 "owner-triangle shallow":([1,2,4,7,12,21,37,65,114],"cubic-Pisot 1.7549 (x^3-2x^2+x-1)"),
}
for name,(seq,const) in ladders.items():
    r=fit(seq); root=charroot(r) if r else None
    print(f"  {name:26s}: {seq[:8]}  -> {const}" + (f"  [rec {r}, root {root:.4f}]" if r and root else ""))
print("\n  THE THREE FAMOUS CONSTANTS (2, golden phi, plastic) = Pascal moved down 0,1,2.")
print("  Higher figurate families + higher slopes birth a whole spectrum of Pisot ratios.")
