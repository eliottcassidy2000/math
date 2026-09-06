#!/usr/bin/env python3
"""Exact finite-dimension hostiles; a two-family control set, not an optimizer census."""
from fractions import Fraction as F
from math import isqrt
import hashlib
import json

gates=0


def check(value,label):
    global gates
    gates+=1
    if not value:
        raise RuntimeError(label)


class I:
    def __init__(self,a,b=None):
        self.a,self.b=F(a),F(a if b is None else b)
        if self.a>self.b:
            raise ValueError('interval order')
    def __add__(self,q):
        q=iv(q)
        return I(self.a+q.a,self.b+q.b)
    __radd__=__add__
    def __neg__(self):
        return I(-self.b,-self.a)
    def __sub__(self,q):
        return self+-iv(q)
    def __rsub__(self,q):
        return iv(q)+-self
    def __mul__(self,q):
        q=iv(q)
        vals=[self.a*q.a,self.a*q.b,self.b*q.a,self.b*q.b]
        return I(min(vals),max(vals))
    __rmul__=__mul__
    def __truediv__(self,q):
        q=iv(q)
        if q.a<=0<=q.b:
            raise ValueError('zero in divisor interval')
        return self*I(1/q.b,1/q.a)
    def __rtruediv__(self,q):
        return iv(q)/self
    def __pow__(self,k):
        result=I(1)
        for _ in range(k):
            result=result*self
        return result
    def contains(self,q):
        return self.a<=q<=self.b
    def display(self):
        scale=10**12
        return [str(F((self.a*scale).__floor__(),scale)),
                str(F((self.b*scale).__ceil__(),scale))]


def iv(q):
    return q if isinstance(q,I) else I(q)


def radical(q):
    q=F(q)
    scale=10**40
    t=isqrt(q.numerator*scale**2//q.denominator)
    result=I(F(t,scale),F(t+1,scale))
    check(result.a**2<=q<=result.b**2,'outward rational square-root certificate')
    return result


U=radical(2)
V=radical(3)
K=4*V/(3*(1+U)*(1+V))
C=-22-F(16,3)*radical(6)+10*U+F(40,3)*V
CSTAR=(13-8*U)/3


def ratio(p3,p4,top_two):
    check((1-iv(p4)).a>0,'positive energy')
    distance=2-U*top_two
    check(distance.a>0,'positive ordered top-two distance')
    J=(5-8*iv(p3)+3*iv(p4))/(3*(1-iv(p4)))
    return (J-CSTAR)/distance


def rational_curve(N,k):
    T=(N-3)*(N-2)//2
    nums=[T*k*k,(N-3)*k,1]+[-k]*(N-3)
    denominator=T*k*k+1
    row=[F(q,denominator) for q in nums]
    check(len(row)==N and all(row),'declared exact nonzero length')
    check(sum(row)==1 and sum(q*q for q in row)==1,'literal exact moment normalization')
    positives=sorted((q for q in row if q>0),reverse=True)
    check(positives[:2]==row[:2],'actual ordered top-two labels')
    r=ratio(sum(q**3 for q in row),sum(q**4 for q in row),sum(positives[:2]))
    return row,r,N*(r-K)


def equal_three(N):
    n=N-3
    a=(1+radical(F(n*(n+2),3)))/(n+3)
    dust=(3*a-1)/n
    check(a.a>0 and dust.a>0,'positive leading triple and negative dust')
    check((3*a-n*dust).contains(1),'first moment interval control')
    check((3*a*a+n*dust*dust).contains(1),'second moment interval control')
    p3=3*a**3-n*dust**3
    p4=3*a**4+n*dust**4
    r=ratio(p3,p4,2*a)
    return r,N*(r-K)


def main():
    check(F(21722,10000)<C.a<C.b<F(21723,10000),'simple certified constant bracket')
    boundary=U-F(2,3)
    check((5*(boundary-K)-C).b<0 and (6*(boundary-K)-C).a>0,
          'one-atom obstruction changes sign between total lengths5 and6')

    witnesses=[]
    for N,first,lo,hi in ((4,16,F(21419,10000),F(21420,10000)),
                         (5,15,F(21649,10000),F(21650,10000))):
        for k in range(1,first+1):
            row,r,penalty=rational_curve(N,k)
            if k<first:
                check((penalty-C).a>0,'earlier integer in the declared curve is not hostile')
            else:
                check((penalty-C).b<0,'actual finite-length counterexample')
                check(lo<penalty.a<penalty.b<hi,'simple rational hostile bracket')
                upper,_=equal_three(N)
                check(r.b<upper.a,'equal-three family is not a finite-N optimizer')
                witnesses.append(dict(N=N,k=k,row=[str(q) for q in row],
                                      R=r.display(),scaled_penalty=penalty.display(),
                                      defect=(C-penalty).display()))

    # Precisely17 dimensions, one fixed rational member and one equal-three member each.
    # This does not search the normalized sphere or any minimizer family.
    controls=[]
    for N in range(4,21):
        row,r,penalty=rational_curve(N,100)
        er,ep=equal_three(N)
        check((penalty-C).b<0 if N<=5 else (penalty-C).a>0,
              'finite coefficient comparison on named rational controls')
        check((ep-C).a>0,'finite coefficient comparison on named equal controls')
        check(r.b<er.a if N<=7 else er.b<r.a,
              'declared two-family comparison, not global optimization')
        controls.append(dict(N=N,rational_R=r.display(),equal_three_R=er.display(),
                             rational_scaled_penalty=penalty.display(),
                             equal_three_scaled_penalty=ep.display()))

    manifest=dict(status='REFUTED all-N coefficient bound and finite-N equal-three optimality',
                  universe='N4..20: two declared controls each; parameter heads k1..16 andk1..15',
                  K3=K.display(),C=C.display(),one_atom_boundary=boundary.display(),
                  witnesses=witnesses,controls=controls,
                  surviving_question='N(R-K3)>=C for every eligible length<=N when N>=6 remains OPEN')
    encoded=json.dumps(manifest,sort_keys=True,separators=(',',':'))
    print('FINITE-EXACT normalized hostiles; asymptotic theorem unchanged')
    print(encoded)
    print('EXPLICIT_GATES',gates)
    print('SEMANTIC_SHA256',hashlib.sha256(encoded.encode()).hexdigest())


if __name__=='__main__':
    main()
