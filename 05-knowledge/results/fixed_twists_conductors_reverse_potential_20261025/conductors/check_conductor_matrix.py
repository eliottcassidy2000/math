# Verify the normalization and conductor formulas for A = R[u1,u2,u3]/(u_i^2 - b_i) over R = Z[1/2] with b = (15, 21, 55):
#   c_z = sqrt-ideal of b^z (p-adic floor(v/2)), A' = sum c_z^{-1} u_z,  conductor component v_p(d_z) = max(0, max_z' (floor(v_p(b^z')/2) - v_p(b^{z AND z'}))).
from fractions import Fraction
from itertools import product
b = {1:15, 2:21, 4:55}            # character bits: z in 0..7, b^z = prod over set bits
primes = [3,5,7,11]
def vp(n,p):
    n=abs(n); c=0
    while n and n%p==0: n//=p; c+=1
    return c
def bz(z):
    r=1
    for bit,val in b.items():
        if z & bit: r*=val
    return r
def cz(z):   # principal generator of the sqrt ideal (R = Z[1/2], only odd primes matter)
    r=1
    for p in primes: r*=p**(vp(bz(z),p)//2)
    return r
# multiplication: u_z u_w = b^{z AND w} u_{z XOR w}
def mul(x,y):
    out={}
    for z,fz in x.items():
        for w,fw in y.items():
            out[z^w]=out.get(z^w,0)+fz*fw*bz(z&w)
    return out
def in_R(q): return Fraction(q).denominator & (Fraction(q).denominator-1) == 0   # denominator a power of 2
def in_A(x): return all(in_R(f) for f in x.values())
gens=[{z:Fraction(1,cz(z))} for z in range(8)]   # generators c_z^{-1} u_z of A'
# (1) closure of A' under multiplication: products of generators lie in A' (coefficient of u_w in c_w^{-1} R)
closed=all(all(in_R(f*cz(w)) for w,f in mul(g,h).items()) for g in gens for h in gens)
print("A' closed under multiplication:", closed)
# (2) integrality: each generator squared lies in R  (and e.g. u_z/(p c_z) is not integral)
print("generators integral (square in R):", all(in_R(mul(g,g).get(0,0)) for g in gens))
print("u_z/(p*c_z) never integral:", all(not in_R(Fraction(bz(z), (p*cz(z))**2)) for z in range(1,8) for p in primes))
# (3) conductor components: f u_z is in (A:A') iff f*g in A for all generators g; compare predicted valuation matrix
ok=True
for z in range(8):
    pred={p: max([0]+[vp(bz(w),p)//2 - vp(bz(z&w),p) for w in range(8)]) for p in primes}
    # minimal f = prod p^{e_p}: test that f u_z is in the conductor and lowering any exponent breaks it
    f=1
    for p in primes: f*=p**pred[p]
    x={z:Fraction(f)}
    inside=all(in_A(mul(x,g)) for g in gens)
    minimal=all((pred[p]==0) or not all(in_A(mul({z:Fraction(f,p)},g)) for g in gens) for p in primes)
    ok &= inside and minimal
    print(f"z={z:03b}: b^z={bz(z)}, c_z={cz(z)}, predicted conductor generator {f}: in conductor {inside}, minimal {minimal}")
print("conductor valuation matrix verified:", ok)
# (4) one-step vs mixed tests: x = u_1 (z=001): single-root tests u_i x in A hold, but the mixed reduced root (u1 u2/3) x ...
x={1:Fraction(1)}
print("u1: single-root products in A:", all(in_A(mul(x,{w:Fraction(1)})) for w in (1,2,4)), "| times u1u2/3 in A:", in_A(mul(x,{3:Fraction(1,3)})))
