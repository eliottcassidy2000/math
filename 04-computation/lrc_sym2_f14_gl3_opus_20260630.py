"""sym^2 f_14 via the EULER PRODUCT (fast): a_p by Legendre-symbol count, L(sym^2,s)=prod_p L_p(s)."""
import math
a1,a2,a3,a4,a6=1,0,1,4,-6
def ap(p):
    if p==2: return -1
    if p==7: return 1
    s=0
    for x in range(p):
        f=(x**3+a2*x*x+a4*x+a6)%p
        disc=(4*f+(a1*x+a3)**2)%p
        s+= 0 if disc==0 else (1 if pow(disc,(p-1)//2,p)==1 else -1)
    return -s
# sieve primes
def primes_up_to(N):
    sieve=bytearray([1])*(N+1);sieve[0]=sieve[1]=0
    for i in range(2,int(N**.5)+1):
        if sieve[i]:
            for j in range(i*i,N+1,i):sieve[j]=0
    return [i for i in range(2,N+1) if sieve[i]]
P=primes_up_to(4000)
apt={p:ap(p) for p in P}
def Lsym2(s):
    prod=1.0
    for p in P:
        A=apt[p]
        if p in (2,7):
            prod*= 1.0/(1-p**(-s))
        else:
            x=p**(-s)
            # 1/[(1-p x)(1-(A^2-2p)x+p^2 x^2)]
            den=(1-p*x)*(1-(A*A-2*p)*x+p*p*x*x)
            prod*=1.0/den
    return prod
print("sym^2 f_14 = the GL(3) degree-3 L-function (the 2nd moment object).")
print("b_p = a_p^2 - p:", {p:apt[p]**2-p for p in [3,5,11,13,17]})
for s in [3.5,3.0,2.7,2.5]:
    print(f"   L(sym^2 f_14, {s}) ~ {Lsym2(s):.5f}   (Euler product, primes<4000)")
print()
floor=114382/332563; bulk=3/math.pi**2; atom=4*math.cos(3*math.pi/7)**2
print(f"   floor: cusp-part {floor-bulk:.5f} | bulk {bulk:.5f} | floor {floor:.5f} | doublet {atom:.5f}")
print()
print("STRUCTURAL (the real payoff, independent of the exact value):")
print("  * the floor is the 2nd MOMENT (pair correlation) of the danger count.")
print("  * L(f x f, s) = zeta(s) . L(sym^2 f, s): Eisenstein bulk = zeta piece, obstruction = sym^2 piece.")
print("  * sym^2 f_14 is degree 3 = GL(3) = SL(3) = the LITTLEWOOD level of the moment hierarchy.")
print("  * so LRC(14)'s obstruction sits one GL-rank ABOVE the curve (GL(2)->GL(3)), exactly where the")
print("    Littlewood/Siegel-Rogers moment thread already placed it (SL(3), n past Littlewood's SL(2)xSL(2)).")
