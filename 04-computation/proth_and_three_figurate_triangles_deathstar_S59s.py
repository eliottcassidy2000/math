#!/usr/bin/env python3
"""
death-star-2026-07-20-S59s (HYP-8165, THM-1355) -- the n*2^x+1 = PROTH unification
and the three figurate triangles (polyhedral / polygonal / power-sum).
"""
from math import comb

print("=== n*2^x+1 = PROTH NUMBER family: the two repo master-constants as axes ===")
def f(x,n): return n*2**x+1
print("  x=1 slice (2n+1) = ODD NUMBERS = observer principle A_k<->n=2k+1, Redei-odd:")
print("   ", [f(1,n) for n in range(1,8)])
print("  n=1 slice (2^x+1) = hypotenuse H=1+2^d / Fermat rungs / CD tower:")
print("   ", [f(x,1) for x in range(1,7)], " (3,5,17,257 = Fermat primes THM-871)")
print("  DIAGONAL x=n (n*2^n+1) = CULLEN NUMBERS:", [f(n,n) for n in range(1,6)])
print("  n=1, x=2^m (2^(2^m)+1) = FERMAT NUMBERS:", [f(2**m,1) for m in range(0,4)])
print("  Proth's theorem: primality of k*2^x+1 (k odd, 2^x>k) is efficiently testable")
print("  => the family the owner wrote is the PROTH numbers; its two axes are the")
print("     repo's observer (2n+1) and hypotenuse (2^x+1). All ODD for x>=1 (Redei sig).")

print("\n=== THE THREE FIGURATE TRIANGLES (columns = a figurate family) ===")
# polyhedral (Pascal simplex): T(n,k)=C(n-1,k); row sums 2^(n-1); col k = simplex
def polyhedral_row(n): return [comb(n-1,k) for k in range(n)]
# polygonal: col j = (j+2)-gonal-ish; opus-S317: row sums = Moser A000127
def moser(n): return 1+comb(n,2)+comb(n,4)
# power-sum (owner's, as given)
owner=[[1],[2,1],[3,3,1],[4,6,5,1],[5,10,14,9,1],[6,15,30,37,17,1],[7,21,55,101,99,33,1]]
print("  POLYHEDRAL (Pascal): row sums = 2^(n-1):", [sum(polyhedral_row(n)) for n in range(1,8)])
print("  POLYGONAL: row sums = Moser A000127     :", [moser(n) for n in range(0,7)], "(1,2,4,8,16,31,57)")
print("  POWER-SUM (owner): row sums             :", [sum(r) for r in owner], "(1,3,7,16,39,106,317)")
print("  TRIANGULAR NUMBERS appear as column 1 in ALL THREE (the 'third perspective'):")
print("    polyhedral col1 = C(n-1,1) = n-1;  polygonal col1;  power-sum col1 = Sum k = triangular")
print("  Owner col1 (triangular):", [owner[r][1] for r in range(1,7)])

print("\n=== THE MOSER-BREAK: owner columns = power sums, exact for j<=2, break at cubes ===")
for j in range(4):
    col=[owner[r][j] for r in range(len(owner)) if j<len(owner[r])]
    ps=[sum(k**j for k in range(1,m+1)) for m in range(1,len(col)+1)]
    status="EXACT" if col==ps else f"BREAKS (col={col} vs Sum k^{j}={ps})"
    print(f"  col j={j}: {status}")
print("  <-> Moser: 1,2,4,8,16 looks like 2^n then breaks (31); hypotenuse 2^(n-2)+1")
print("     is exact single-flip but BREAKS at n=7 for full tiling class (THM/tiling-class).")
print("     THREE independent 'looks-simple-then-breaks' at the SAME lesson.")

print("\n=== FIBONACCI: shallow diagonals ===")
def shallow(tri,d):
    return sum(tri[d-1-k][k] for k in range(d) if 0<=d-1-k<len(tri) and k<len(tri[d-1-k]))
print("  Pascal shallow diagonals -> FIBONACCI (x^2=x+1, phi=1.618)")
owner_shallow=[shallow(owner,d) for d in range(1,8)]
print("  owner-triangle shallow diagonals:", owner_shallow, " (cubic-Fibonacci)")
print("  recurrence a(n)=2a(n-1)-a(n-2)+a(n-3); growth root x^3-2x^2+x-1 ~ 1.7549")
print("  => the power-sum triangle's 'Fibonacci' is a CUBIC Pisot analog of the golden ratio.")

print("\n=== HARMONIC (series 1) = the STAIRCASE electrical resistance R_{n-2} (HYP-6865) ===")
from fractions import Fraction as Fr
H=[sum(Fr(1,k) for k in range(1,m+1)) for m in range(1,6)]
print("  H_n:", [str(h) for h in H], " = R_{n-2} pole resistance of staircase Delta_{n-2}; gamma = growth")
print("  series 2 (2-power harmonic Sum (2^k-1)/k):", [str(sum(Fr(2**j-1,j) for j in range(1,m+1))) for m in range(1,6)])
print("  the two series = log-shadow (harmonic~ln n) and exp-shadow (~2^n/n): the two axes again.")
