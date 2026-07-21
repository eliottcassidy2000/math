#!/usr/bin/env python3
"""
klein-2026-07-20-S381 -- IS THE GMC(2) HOLONOMIC ORDER OF E[P^m] BOUNDED BY THE SPAN, uniform
in radial degree?  (HYP-8540 = the radial analog of THM-1710's toral 'depth = span'.)

E[P^m] is P-recursive in m (opus/kp THM-1740).  If its recurrence ORDER is bounded by the SPAN
alone -- uniformly over the radial degree d -- then with a non-degenerate leading coefficient
E[P^1]=..=E[P^D]=0 forces all moments to vanish, so GMC(2) holds for ALL bounded-span P at every
degree.  This measures that order by finding the minimal P-finite recurrence of the exact
integer moment sequence, for genuine (charge-radius-LOCKED) two-sided P.
"""
from fractions import Fraction as Fr
from math import factorial
import random

def moments_numeric(charges, dmax, coeff, M):
    """E[P^m], m=1..M, for the locked P with numeric rational coeff[(q,k)].
       monomials: charge q, radial depth k -> Z^a Zb^b, a=(|q|+2k+q)/2, b=(|q|+2k-q)/2."""
    # represent P as dict {(a,b): coeff}
    P = {}
    for q in charges:
        for k in range(dmax+1):
            h = abs(q)+2*k; a=(h+q)//2; b=(h-q)//2
            P[(a,b)] = P.get((a,b),Fr(0)) + coeff[(q,k)]
    # power via dict convolution, then E = sum over (a,a) of coeff*a!
    def mul(A,B):
        C={}
        for (a1,b1),c1 in A.items():
            for (a2,b2),c2 in B.items():
                key=(a1+a2,b1+b2); C[key]=C.get(key,Fr(0))+c1*c2
        return C
    out=[]; cur={(0,0):Fr(1)}
    for m in range(1,M+1):
        cur=mul(cur,P)
        e=Fr(0)
        for (a,b),c in cur.items():
            if a==b: e+=c*factorial(a)
        out.append(e)
    return out

def min_recurrence_order(seq, max_ord, max_deg):
    """smallest order D such that some P-finite recurrence sum_i c_i(m) a_{m+i}=0 holds,
       c_i poly in m of degree <= max_deg.  Returns D or None."""
    M=len(seq)
    for D in range(1, max_ord+1):
        for e in range(0, max_deg+1):
            nun=(D+1)*(e+1)
            rows=[]
            for m in range(1, M-D+1):          # equation at shift m
                row=[]
                for i in range(D+1):
                    for j in range(e+1):
                        row.append(Fr(m)**j * seq[m-1+i])
                rows.append(row)
            if len(rows) < nun+2: continue     # need enough equations to avoid spurious
            # null space over Q by gaussian elimination
            if has_nontrivial_null(rows, nun):
                return D
    return None

def has_nontrivial_null(rows, ncol):
    # rank via exact gaussian elimination; nontrivial null iff rank < ncol
    M=[r[:] for r in rows]; nr=len(M); rank=0; col=0
    while col<ncol and rank<nr:
        piv=next((r for r in range(rank,nr) if M[r][col]!=0), None)
        if piv is None: col+=1; continue
        M[rank],M[piv]=M[piv],M[rank]
        inv=1/M[rank][col]
        M[rank]=[x*inv for x in M[rank]]
        for r in range(nr):
            if r!=rank and M[r][col]!=0:
                f=M[r][col]; M[r]=[a-f*b for a,b in zip(M[r],M[rank])]
        rank+=1; col+=1
    return rank<ncol

print("="*92)
print("HOLONOMIC ORDER of E[P^m] vs (charge span, radial degree d) -- for random locked two-sided P")
print("="*92)
print(f"{'charges':>18} {'span':>5} {'d':>3} {'M':>4} {'holonomic order (2 random P)':>30}")
random.seed(3)
tests=[([-1,0,1],0),([-1,0,1],1),([-1,0,1],2),([-1,0,1],3),
       ([-1,1],0),([-1,1],1),([-1,1],2),
       ([-2,-1,1,2],0),([-2,-1,1,2],1),
       ([-2,-1,0,1,2],0),([-2,-1,0,1,2],1),
       ([-2,0,2],0),([-3,-1,1,3],0)]
for charges,d in tests:
    span=max(charges)-min(charges)
    orders=[]
    for trial in range(2):
        coeff={(q,k):Fr(random.randint(-5,5) or 1, random.randint(1,3)) for q in charges for k in range(d+1)}
        M=40
        seq=moments_numeric(charges,d,coeff,M)
        D=min_recurrence_order(seq, max_ord=span+4, max_deg=span+d+4)
        orders.append(D)
    print(f"{str(charges):>18} {span:>5} {d:>3} {M:>4} {str(orders):>30}")
print("""
 READING: THM-1710 gives toral (d=0) order = span.  If the order stays = span as d grows
 (constant down each fixed-charges block), the HOLONOMIC ORDER IS SPAN-UNIFORM -> HYP-8540 holds
 -> GMC(2) for ALL bounded-span P at every radial degree.  If it grows with d, the radial layer
 needs strictly more than the toral depth and HYP-8540 (span-only bound) is FALSE.
""")
