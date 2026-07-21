#!/usr/bin/env python3
"""klein-2026-07-20-S383 -- GROWTH LAW of the GMC(2) holonomic order vs radial degree, mod p
(fast; order is prime-independent for a generic prime). Load-bearing follow-up to THM-1770."""
import random
p = 2_000_000_011                       # prime > any factorial argument used (max radius ~ 50*8)
def moments(charges, dmax, coeff, M):
    P = {}
    for q in charges:
        for k in range(dmax+1):
            h=abs(q)+2*k; a=(h+q)//2; b=(h-q)//2
            P[(a,b)]=(P.get((a,b),0)+coeff[(q,k)])%p
    # factorials mod p
    fmax=0
    # build power dict
    def mul(A,B):
        C={}
        for (a1,b1),c1 in A.items():
            for (a2,b2),c2 in B.items():
                key=(a1+a2,b1+b2); C[key]=(C.get(key,0)+c1*c2)%p
        return C
    out=[]; cur={(0,0):1}
    fact=[1]*(1)
    for m in range(1,M+1):
        cur=mul(cur,P)
        # extend factorial table
        mx=max((a for (a,b) in cur if a==b), default=0)
        while len(fact)<=mx: fact.append(fact[-1]*len(fact)%p)
        e=0
        for (a,b),c in cur.items():
            if a==b: e=(e+c*fact[a])%p
        out.append(e)
    return out
def rref_null(rows, ncol):
    M=[r[:] for r in rows]; nr=len(M); piv={}; r=0
    for c in range(ncol):
        pr=next((i for i in range(r,nr) if M[i][c]!=0), None)
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]; inv=pow(M[r][c],p-2,p)
        M[r]=[(x*inv)%p for x in M[r]]
        for i in range(nr):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],M[r])]
        piv[c]=r; r+=1
    free=[c for c in range(ncol) if c not in piv]
    if not free: return None
    f=free[0]; x=[0]*ncol; x[f]=1
    for c,rr in piv.items(): x[c]=(-M[rr][f])%p
    return x
def order(seq, max_ord, max_deg):
    M=len(seq)
    for D in range(1,max_ord+1):
        for e in range(0,max_deg+1):
            nun=(D+1)*(e+1); K=nun+D+3
            if K+D+6>M: continue
            rows=[[ (pow(m,j,p)*seq[m-1+i])%p for i in range(D+1) for j in range(e+1)] for m in range(1,K+1)]
            x=rref_null(rows,nun)
            if x is None: continue
            ok=all(sum(x[i*(e+1)+j]*pow(m,j,p)*seq[m-1+i] for i in range(D+1) for j in range(e+1))%p==0
                   for m in range(K+1,M-D+1))
            if ok: return D
    return None
print("="*78); print("HOLONOMIC ORDER of E[P^m] vs (span, radial degree d) -- mod p, cross-validated, M=55"); print("="*78)
random.seed(11); M=55
print(f"{'charges':>16} {'span':>5} {'d':>3} {'order (3 generic P)':>22}")
for charges in ([-1,0,1],[-1,1],[-2,-1,0,1,2]):
    for d in range(0,4):
        ords=[]
        for _ in range(3):
            coeff={(q,k):random.randint(1,p-1) for q in charges for k in range(d+1)}
            seq=moments(charges,d,coeff,M)
            ords.append(order(seq, max_ord=4*(d+1)+2, max_deg=3*(d+1)+3))
        print(f"{str(charges):>16} {max(charges)-min(charges):>5} {d:>3} {str(ords):>22}")
print("""
 THM-1710/1770: order = span at d=0.  A clean law in d is the 'resonance complexity' the bridge
 must dominate -- load-bearing.  (3 generic P per row; agreement across them = reliable order.)
""")
