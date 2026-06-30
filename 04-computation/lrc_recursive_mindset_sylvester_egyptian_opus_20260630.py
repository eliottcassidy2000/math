"""
Push the a=x/2 (descend), b=x+1 (observer) recursion everywhere.
KEY: iterating Phi_6(n)=n(n-1)+1 = b(2T_{n-1}) IS Sylvester's sequence (2,3,7,43,1807,...) = the greedy
Egyptian fraction of the observer's 1 = 1/2+1/3+1/7+1/43+... The apex primes 3,7 are its early terms.
Also: the doublet {0,1}={0,b(0)}=origin+observer (apex attractor); Cayley-Dickson n=2^k+1=b(2^k); fixed pt 1.
"""
from fractions import Fraction
Phi6=lambda n: n*n-n+1
b=lambda x:x+1; a=lambda x:x//2
# (1) the Phi_6-iteration = Sylvester
syl=[2]
for _ in range(5): syl.append(Phi6(syl[-1]))
print(f"Phi_6-iteration n->n(n-1)+1 from 2 = SYLVESTER's sequence: {syl}")
print(f"  = the apex primes 3,7 (early terms) + 43,1807,...; Phi_6(2)=3, Phi_6(3)=7, Phi_6(7)=43.")
# (2) Egyptian fraction of the observer's 1 (Sylvester greedy) + telescoping 1/(a_k-1)-1/(a_{k+1}-1)=1/a_k
S=sum(Fraction(1,t) for t in syl[:5])
print(f"\n  observer's 1 = sum 1/a_k (Sylvester): 1/2+1/3+1/7+1/43+1/1807+... = {float(S):.8f} -> 1")
print(f"  telescoping: 1/(a_k - 1) - 1/(a_{{k+1}}-1) = 1/a_k  (since a_{{k+1}}-1 = a_k(a_k-1)):")
for i in range(4):
    lhs=Fraction(1,syl[i]-1)-Fraction(1,syl[i+1]-1); print(f"    1/{syl[i]-1} - 1/{syl[i+1]-1} = {lhs} = 1/{syl[i]}  {lhs==Fraction(1,syl[i])}")
# (3) the doublet {0,1} = {origin, observer}; apex attractor
print(f"\n  DOUBLET {{0,1}} = {{0, b(0)}} = origin + observer = the descent attractor / apex minimal core.")
print(f"     the apex gap 4cos^2(3pi/7) lives on this {{0,b(0)}} pair; the +1 IS the observer in the core.")
# (4) Cayley-Dickson tower n=2^k+1=b(2^k); and apex primes 3,7 = b(2)=b(2T_1), b(6)=b(2T_2)
print(f"\n  CAYLEY-DICKSON tower R,C,H,O,S: n = 2^k+1 = b(2^k): {[2**k+1 for k in range(5)]}  (=b on the doubling).")
print(f"     apex primes: 3=b(2T_1)=Phi_6(2), 7=b(2T_2)=Phi_6(3) -- the observer +1 on doubled triangles.")
# (5) fixed points: everything flows to 1 (the observer)
print(f"\n  FIXED POINT 1 (the observer's irreducible baseline):")
print(f"     a.b (=g on odds) (x+1)/2 fixes 1; Phi_6(1)=1; the descent g^k(x)->1. ALL flows converge to 1.")
# (6) the two columns, both b(+1)-structured
print(f"\n  TWO COLUMNS, both built by b(+1):")
print(f"     MEASURE/apex Q(-7): doublet {{0,1}}={{0,b(0)}} (minimal core)  -> gap 4cos^2(3pi/7)")
print(f"     EXIST/cover Q(-3):  Phi_6=2T_{{n-1}}+1=b(2T) (max construction) -> M=n/Phi_6")
print(f"     => b(+1) at BOTH ends: the observer on the minimal core (doublet) and on the doubled triangle (cover).")
