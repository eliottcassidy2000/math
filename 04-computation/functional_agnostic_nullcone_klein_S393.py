#!/usr/bin/env python3
"""
klein-2026-07-20-S393 -- THE FUNCTIONAL-AGNOSTIC BOTH-SIGNS SINGLE-CHARACTER NULLCONE
NON-VANISHING, and its cyclotomic / coprime-pair-return-time content.  The rigorous common
parent of THM-1830.

Owner: "prove the abstract 'both-signs single-character nullcone non-vanishing' by the
functional-agnostic method; think cyclotomic."

THEOREM (single-character = two charges).  Let Lambda(u) = a u^p + b u^{-q}, a,b != 0, p,q>0
(BOTH-SIGNS, one positive charge +p and one negative charge -q).  Let g = gcd(p,q),
m0 = (p+q)/g.  Then the CHARGE-0 coefficient of Lambda^m is

    [u^0] Lambda^m = C(m, qm/(p+q)) a^{qm/(p+q)} b^{pm/(p+q)}   if  m0 | m,   else 0,

a SINGLE term (the only (k,m-k) split with p k = q(m-k)), hence UNCANCELLABLE and NONZERO at
every return time m in m0*Z_{>0}.  THEREFORE, for ANY charge-graded functional F (one that
depends on a monomial only through its charge, F(u^0) != 0), F(Lambda^m) = F(u^0) * [u^0]Lambda^m
!= 0 at those m -- the non-vanishing is a fact about the CHARGE LATTICE, independent of which
functional (CT/TNC, Gaussian/GMC, Bonferroni/LRC) evaluates it.  m0 = (p+q)/gcd(p,q) is the
COPRIME-PAIR FIRST-RETURN TIME (mac-mini THM-1745) and the condition p k = q(m-k) is CYCLOTOMIC
(balance on the charge lattice = the (p+q)-th roots of unity).
"""
import itertools
from math import comb, gcd, factorial
from fractions import Fraction as Fr

def charge0_binomial(a,b,p,q,m):
    """[u^0] of (a u^p + b u^{-q})^m -- brute, to verify the single-term formula"""
    tot=0; terms=0
    for k in range(m+1):
        if p*k - q*(m-k) == 0:
            tot += comb(m,k)*a**k*b**(m-k); terms+=1
    return tot, terms

print("="*84)
print("(1) BINOMIAL: the charge-0 coefficient is a SINGLE term at the return times")
print("="*84)
print(f"{'(p,q)':>8} {'m0=(p+q)/g':>12} {'m':>4} {'[u^0]Lam^m (a=2,b=3)':>22} {'#terms':>7}")
for (p,q) in [(1,1),(1,2),(2,3),(3,5),(2,2)]:
    g=gcd(p,q); m0=(p+q)//g
    for m in (m0, 2*m0, m0+1):
        val,terms=charge0_binomial(2,3,p,q,m)
        note = "  <- return time, SINGLE nonzero term" if (m%m0==0 and val!=0) else ("  (not a return time)" if m%m0 else "")
        print(f"{str((p,q)):>8} {m0:>12} {m:>4} {str(val):>22} {terms:>7}{note}")
print("""
 The charge-0 term is present ONLY when m0=(p+q)/gcd(p,q) divides m, and then it is EXACTLY ONE
 term -- uncancellable.  m0 is the coprime-pair first-return time.
""")

print("="*84)
print("(2) FUNCTIONAL-AGNOSTIC: three DIFFERENT charge-graded functionals, same non-vanishing")
print("="*84)
print(" Lambda = 2 u^2 + 3 u^{-3}  (p=2,q=3,m0=5).  Evaluate F(Lambda^5) for:")
print("   F_CT   : the constant-term (toral/TNC) functional, F(u^j)=delta_{j0}          [angular]")
print("   F_gau  : charge-graded Gaussian-type, F(u^j)=delta_{j0}*w0 with w0 = 7 (a radial weight)")
print("   F_bonf : charge-graded 'Bonferroni-type', F(u^j)=delta_{j0}*(-1)^?*W, W = 5")
p,q=2,3; m0=5; a,b=2,3
c0,_=charge0_binomial(a,b,p,q,m0)
print(f"   [u^0]Lambda^5 = {c0}  (the charge fact; = C(5,3) a^3 b^2 = {comb(5,3)*a**3*b**2})")
for name,w0 in [("F_CT",1),("F_gau",7),("F_bonf",5)]:
    print(f"     {name}(Lambda^5) = w0 * [u^0] = {w0}*{c0} = {w0*c0}  (nonzero for ANY w0 != 0)")
print("""
 The non-vanishing is IDENTICAL up to the scalar F(u^0): the functional supplies a nonzero
 weight, the CHARGE LATTICE supplies the non-vanishing.  This is the functional-agnostic method
 -- it proves the both-signs single-character nullcone non-vanishing simultaneously for TNC
 (F=CT), GMC (F = radial-weighted), and LRC (F = Bonferroni-weighted), because all three are
 charge-graded and faithful on charge 0.
""")

print("="*84)
print("(3) GENERAL BOTH-SIGNS: the extreme charges give a return-time term -- unique iff binomial")
print("="*84)
print(" For Lambda with charges in [-Q,+P], the extreme pairing (+P with -Q) yields a charge-0")
print(" term at m0=(P+Q)/gcd(P,Q); it is UNIQUE (uncancellable) exactly in the binomial case.")
print(" With intermediate charges present it may COLLIDE with other charge-0 terms -- and THAT")
print(" collision is exactly the GMC-hard many-charge cancellation (THM-1810 bosonic side).")
# demonstrate: {-1,0,1} has extreme pair (+1,-1), m0=2, but the charge-0 part of Lam^2 also has
# the (0,0) self-term -> two contributions, can cancel.
def charge0_general(coeffs, m):
    """coeffs: dict charge->coeff.  [u^0] of (sum c_q u^q)^m, brute."""
    from collections import defaultdict
    cur={0:1}
    for _ in range(m):
        nxt=defaultdict(int)
        for e,c in cur.items():
            for qy,cq in coeffs.items(): nxt[e+qy]+=c*cq
        cur=nxt
    return cur.get(0,0)
print("   Lambda = 1*u^1 + B*u^0 + 1*u^{-1}, m=2: [u^0] = 2*1*1 + B^2 = 2 + B^2")
print(f"     B=0: {charge0_general({1:1,0:0,-1:1},2)}   B such that 2+B^2=0 => B=i sqrt2 (the THM-1770 witness!)")
print("   => the EXTREME-pair term (2) is CANCELLABLE by the charge-0 self-term (B^2) once a")
print("      middle charge exists.  Binomial (no middle charge) is exactly the uncancellable case.")
print("""
 SO: the functional-agnostic method CLOSES the both-signs SINGLE-CHARACTER (binomial) nullcone
 for every functional at once, and PINPOINTS the residual as many-charge cancellation -- the
 same wall THM-1810 reads as 'bosonic/permanent, no sign to cancel'.  The clean cases are
 cyclotomic (return times); the hard cases are where several cyclotomic returns collide.
""")
