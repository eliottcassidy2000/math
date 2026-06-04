from fractions import Fraction as F
from itertools import combinations
from math import comb

# ===== q-ary Krawtchouk: K^(q)_k(n,x) = sum_j (-1)^j (q-1)^{k-j} C(x,j) C(n-x,k-j) =====
def Kq(q,k,n,x): return sum((-1)**j*(q-1)**(k-j)*comb(x,j)*comb(n-x,k-j) for j in range(k+1))
print("=== q-ary Krawtchouk (q=2 -> binary, q=3 -> ternary) ===")
print("  reduces to binary at q=2:", all(Kq(2,k,n,x)==sum((-1)**j*comb(x,j)*comb(n-x,k-j) for j in range(k+1)) for n in range(6) for x in range(n+1) for k in range(n+1)))
print("  Kq_at_zero: K^(q)_k(n,0) = (q-1)^k C(n,k):", all(Kq(q,k,n,0)==(q-1)**k*comb(n,k) for q in[2,3,4] for n in range(7) for k in range(n+1)))
print("  Kq_zero_index: K^(q)_0(n,x)=1:", all(Kq(q,0,n,x)==1 for q in[2,3,4] for n in range(6) for x in range(n+1)))
print("  genfun: sum_k Kq(q,k,n,x) z^k = (1-z)^x (1+(q-1)z)^{n-x}  [q=3 example n=4]:")
for x in range(5):
    coeffs=[Kq(3,k,4,x) for k in range(5)]
    print(f"    x={x}: K^(3)_k = {coeffs}   (= coeffs of (1-z)^{x}(1+2z)^{4-x})")

# ===== n=14: the 3-structure. 2 has ORDER 3 mod 7 -> doubling orbits are 3-cycles =====
print("\n=== n=14 = 2*7: the built-in 3-structure (2 has order 3 mod 7) ===")
def orbit(gen,mod,start):
    o=[start]; x=(gen*start)%mod
    while x!=start: o.append(x); x=(gen*x)%mod
    return o
print("  doubling <2> orbits on (Z/7)*:", orbit(2,7,1), orbit(2,7,3), " (two 3-cycles; ord_7(2)=3)")
print("  antipodal <-1> orbits (order 2):", [(a,(7-a)%7) for a in [1,2,3]], " -> 2-structure (sigma) sits INSIDE the 3-cycles")
print("  so n=14 carries BOTH: 2-structure (sigma pairs) and 3-structure (<2> triples) = the 2-adic/3-adic seam")

# ===== ternary/signed depth refinement: each runner is +(in (0,delta)), -(in (-delta,0)), or 0(safe) =====
def signed_states(V,delta,t):
    sp=sm=0
    for v in V:
        x=(v*t)%1
        d=min(x,1-x)
        if d<delta:
            if x<delta: sp+=1     # near origin from above: + side
            else: sm+=1            # near 1 (=-): - side
    return sp,sm
delta=F(1,14); V=[1,2,3,4,5,6,7,8,9,10,11,13,14]
print("\n=== signed (ternary) depth refinement at n=14 unit clocks: cross = (n+,n-)=(1,1) ===")
for j in [1,3,5,9,11,13]:
    sp,sm=signed_states(V,delta,F(j,14))
    print(f"  clock {j}/14: (n+,n-)=({sp},{sm})  {'<-- CROSS (1,1): one + one -, sigma swaps them' if (sp,sm)==(1,1) else ''}")
print("  binary depth = n+ + n- loses the +/- split; the TERNARY enumerator keeps it (sigma-symmetric: n+ <-> n-)")
print("  lonely <=> (n+,n-)=(0,0); the ternary LP has MORE constraints (n+,n- separately, sigma-symmetry) than the binary LP")

# ===== the 2->3 dictionary =====
print("\n=== the 2 -> 3 dictionary (more concepts like the diagonal certificate, but '3') ===")
for a,b in [
 ("binary Krawtchouk K_k (q=2)","ternary Krawtchouk K^(3)_k (factor 2^{k-j})"),
 ("depth: forbidden / safe (2 states)","signed depth: + / 0 / - (3 states)"),
 ("Helly number 2 (intervals on a LINE)","Helly number 3 (arcs on a CIRCLE) [HYP-2200]"),
 ("pairwise overlap S_2","triple overlap S_3 = 3-term coincidence = a+b=c [HYP-2195]"),
 ("2-adic doubling (div 2)","3-adic / Collatz times-3 [HYP-2175]"),
 ("antipodal involution <-1>, order 2","doubling <2> mod 7, order 3 (the 3-cycles)"),
 ("apex/2-block at n=2q","3-block at n=3q (ternary apex)"),
 ("sigma-pair {v,-v} (the cross)","<2>-triple {v,2v,4v} (the 3-cycle)"),
]:
    print(f"  2: {a:42}  ->  3: {b}")
