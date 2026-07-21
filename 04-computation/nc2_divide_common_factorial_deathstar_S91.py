#!/usr/bin/env python3
"""nc2_divide_common_factorial_deathstar_S91.py (HYP-8795)
IDEA (owner): divide E[P^m] by the full common factorial. Develop the normalization
that boxeph (det/prod a_i! = Vandermonde) and S90 (E[P^m]/m! -> central trinomial)
used, for the general three-weight P=Z^p a(s)+b(s)+Zbar^q c(s). Identify the dominant
common factorial, divide, study the normalized channel sum: does it reveal a clean
tournament (Vandermonde/discriminant) core + sharpen noncancellation?

AUDIT: the interpretation formerly printed at the end is false; see MISTAKE-215.
This script is retained as a counterexample-producing audit of that interpretation.
"""
from math import factorial as f
from fractions import Fraction as Fr

def Lval(coeffs):  # L(sum c_k s^k) = sum c_k k!
    return sum(v*f(k) for k,v in coeffs.items())
def pm(p,q):
    r={}
    for k,v in p.items():
        for k2,v2 in q.items(): r[k+k2]=r.get(k+k2,0)+v*v2
    return {k:v for k,v in r.items() if v}
def pw(p,n):
    r={0:Fr(1)}
    for _ in range(n): r=pm(r,p)
    return r

def EPm_general(p,q,a,b,c,m):
    # P=Z^p a(s)+b(s)+Zbar^q c(s); charge 0: p*j = q*k ; g=gcd; j=qt/g,k=pt/g
    from math import gcd
    g=gcd(p,q); tot=Fr(0)
    t=0
    while True:
        j=q*t//g; k=p*t//g
        if j+k>m: break
        # multinomial m!/(j! k! (m-j-k)!) ; radial: Z^{pj}W^{qk}=s^{pj} (pj=qk); times a^j c^k b^{m-j-k}
        pj=p*j
        multi=Fr(f(m),f(j)*f(k)*f(m-j-k))
        rad=pm(pm(pw(a,j),pw(c,k)),pm(pw(b,m-j-k),{pj:Fr(1)}))
        tot+=multi*Lval(rad)
        t+=1
    return tot

print("Normalizing E[P^m] by the dominant common factorial: what is revealed?\n")
# single-straddle p=q=1, constant a=c=1, b=1+s (a degree-gap-ish): dominant factorial
print("[A] p=q=1, a=c=1, b=1+s:  E[P^m] and E[P^m]/(m!) :")
a={0:Fr(1)}; c={0:Fr(1)}; b={0:Fr(1),1:Fr(1)}
for m in range(1,9):
    E=EPm_general(1,1,a,b,c,m)
    print(f"  m={m}: E[P^m]={E}, /m! = {Fr(E,f(m))}  (float {float(E)/f(m):.4f})")
print("  -> AUDIT: these displayed values grow, so they do NOT support boundedness.\n")

# ASYMMETRIC p=2,q=1: the (pA_0)! idea -- dominant factorial from the busier (Z^p) charge
print("[B] p=2,q=1, a=1,c=1,b=1+s (asymmetric top): E[P^m] / (dominant factorial):")
for m in range(1,10):
    E=EPm_general(2,1,{0:Fr(1)},{0:Fr(1),1:Fr(1)},{0:Fr(1)},m)
    # dominant channel j=m//3 uses Z^{2j}, top factorial ~ (2*(m//3))! -- the 'p*A_0' with A_0=#a-uses
    jmax=m//3; domfac=f(2*jmax) if jmax else 1
    print(f"  m={m}: E={E}, /(2*floor(m/3))!={Fr(E,domfac)} (float {float(E)/domfac:.3f}); jmax={jmax}")

wall_checks = 0
for even_p in (2, 4, 6):
    for aa, beta, cc in ((1, 1, 1), (2, -3, 5), (-1, 2, 3)):
        for m in range(1, 9):
            direct = EPm_general(
                even_p,
                even_p,
                {0: Fr(aa)},
                {even_p // 2: Fr(beta)},
                {0: Fr(cc)},
                m,
            )
            residue = sum(
                Fr(f(m), f(i) * f(i) * f(m - 2 * i))
                * Fr(aa * cc) ** i
                * Fr(beta) ** (m - 2 * i)
                for i in range(m // 2 + 1)
            )
            assert direct == f(even_p * m // 2) * residue
            wall_checks += 1
print(f"\n[C] THM-2040 corrected symmetric-monomial identities checked: {wall_checks} PASS")
print("""
AUDIT (MISTAKE-215):
 * The output above does not establish a common-factorial factorization or a bounded
   normalized sequence.
 * A scalar moment sum is not the moment-matrix determinant in THM-2033, so dividing
   by one scalar factorial does not generally leave a Vandermonde.
 * The valid general normalization is prime-local: THM-2022 chooses one amplified
   moment p*m0 and divides by the lowest-face factorial (p*A0)!; modulo p, carries and
   strict face height kill all other layers and the residue is Q^p.
 * The exact central-trinomial factorization survives only on the symmetric monomial
   wall where every balanced channel has the same radial degree.""")
