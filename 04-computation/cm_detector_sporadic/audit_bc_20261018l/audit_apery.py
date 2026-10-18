#!/usr/bin/env python3
"""Independent all-coefficient check of reflected Apéry carry descent."""
from math import comb
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_identities import primes,trim,add,scale,mul
from audit_zagier import sequence
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ordinary_detector_and_mechanism import mechanism

def apery(N):
    u=[1,5]
    for n in range(2,N+1):
        num=(2*n-1)*(17*n*(n-1)+5)*u[n-1]-(n-1)**3*u[n-2]
        v,r=divmod(num,n**3);assert r==0
        u.append(v)
    return u

def homogeneous(f,A,B,m,p):
    # B^m f(A/B), using Horner and precomputed denominator powers.
    c=list(f)+[0]*(m+1-len(f));q=[c[m]%p];bp=[1]
    for k in range(m-1,-1,-1):
        bp=mul(bp,B,p);q=add(mul(q,A,p),scale(bp,c[k],p),p)
    return q

def packet(u,p):
    f=trim([u[k]%p for k in range(p)])
    ds=[u[p+k]-u[1]*u[k] for k in range(p)]
    assert all(v%p==0 for v in ds)
    v=trim([v//p%p for v in ds])
    g=scale(v,pow(u[1],-1,p),p)
    return f,add(g,scale(f,-g[0],p),p),v

def one(p,uf,ua):
    H,G,Vf=packet(uf,p);F,Ga,Va=packet(ua,p);eps=1 if p%6==1 else p-1
    K=[comb(p-1,k)%p for k in range(p)]
    q=[1,1];sn=[pow(8,-1,p),p-1]
    R_H=homogeneous(H,sn,q,p-1,p)
    R_G=homogeneous(G,sn,q,p-1,p)
    assert H==scale(R_H,eps,p),(p,'H reflection')
    evG=sum(G[k]*pow(8,-k,p) for k in range(len(G)))%p
    B=add(add(G,scale(R_G,eps,p),p),scale(H,-eps*evG,p),p)
    assert B[0]==0
    assert B==scale(homogeneous(B,sn,q,p-1,p),eps,p),(p,'B reflection')
    tn=[0,1,p-8]
    lhsF=homogeneous(F,tn,q,p-1,p)
    lhsG=homogeneous(Ga,tn,q,p-1,p)
    assert lhsF==mul(H,H,p),(p,'period descent')
    assert lhsG==mul(H,B,p),(p,'actual carry descent')
    # Constant normalization is tested, not needed as an input to the proof.
    assert Vf[0]==0 and Va[0]==0
    mech=mechanism(ua,p,True);assert mech['rad|V']
    return {'p':p,'degree_H':len(H)-1,'degree_B':len(B)-1,'degree_F_A':len(F)-1,'degree_G_A':len(Ga)-1,'cleared_period_identity':True,'cleared_carry_identity':True,'reflection_equivariance':True,'mechanism':mech}

if __name__=='__main__':
    bound=int(sys.argv[1]) if len(sys.argv)>1 else 397
    uf=sequence(7,2,-8,2*bound);ua=apery(2*bound)
    assert all(ua[n]==sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1)) for n in range(101))
    out={'scope':'Every prime7<=p<=bound; independent integer recurrence and exact homogeneous substitution; existing mechanism(). No evaluation-only check.','bound':bound,'direct_binomial_control':'All n=0..100 for Apery; Franel separately in audit_zagier.','rows':[one(p,uf,ua) for p in primes(bound) if p>=7]}
    # Retaining the false reciprocal factor fails as a polynomial identity.
    p=7;H,_,_=packet(uf,p);F,_,_=packet(ua,p);K=[comb(p-1,k)%p for k in range(p)]
    cleared=homogeneous(F,[0,1,p-8],[1,1],p-1,p)
    wrong=mul(mul(K,K,p),mul(H,H,p),p)
    assert cleared!=wrong
    out['hostile_controls']={'reciprocal_gauge_error_detected_p7':True}
    print(json.dumps(out,indent=2))
