#!/usr/bin/env python3
"""Independent exact controls for the integral full-sedenion shell corollary.
Prints a deterministic receipt; writes no files. Finite checks are not the proof.
"""
import json
from math import isqrt

GATES = 0
def gate(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)
def add(x,y): return tuple(a+b for a,b in zip(x,y))
def neg(x): return tuple(-a for a in x)
def sub(x,y): return add(x,neg(y))
def scale(k,x): return tuple(k*a for a in x)
def conj(x): return (x[0],)+tuple(-a for a in x[1:])
def dot(x,y): return sum(a*b for a,b in zip(x,y))
def mul(x,y):
    if len(x)==1: return (x[0]*y[0],)
    h=len(x)//2
    a,b=x[:h],x[h:]; c,d=y[:h],y[h:]
    return sub(mul(a,c),mul(conj(d),b))+add(mul(d,a),mul(b,conj(c)))
def four(m):
    for p in range(isqrt(m)+1):
        for q in range(isqrt(m-p*p)+1):
            for r in range(isqrt(m-p*p-q*q)+1):
                s=isqrt(m-p*p-q*q-r*r)
                if p*p+q*q+r*r+s*s==m: return p,q,r,s
    raise RuntimeError('four-square search failed')
def test_shell(m, coefficients):
    p,q,r,s=coefficients
    a=(0,p,q,r,s,0,0,0); b=(0,-q,p,-s,r,0,0,0)
    gate(dot(a,a)==m and dot(b,b)==m,'equal norms')
    gate(dot(a,b)==0,'orthogonality')
    x=a+b
    gate(dot(x,x)==2*m,'full norm shell')
    # Integral projection onto span(1,a,b,ab)^perp, clearing m^2.
    ab=mul(a,b)
    c=None
    for j in range(1,8):
        e=tuple(int(i==j) for i in range(8))
        candidate=sub(sub(sub(scale(m*m,e),scale(m*a[j],a)),scale(m*b[j],b)),scale(ab[j],ab))
        if any(candidate): c=candidate; break
    gate(c is not None and c[0]==dot(c,a)==dot(c,b)==dot(c,ab)==0,'nonzero perpendicular integer c')
    y=scale(m,c)+mul(b,mul(c,a))
    gate(any(y) and not any(mul(x,y)),'explicit integral zero product')
    # Stable higher-stage embedding, preserving this zero product.
    gate(not any(mul(x+(0,)*16,y+(0,)*16)),'embedded zero product at32')
    return x,y

def parity_product(x,y):
    out=[0]*len(x)
    for i,a in enumerate(x):
        for j,b in enumerate(y): out[i^j]^=(a&1)&(b&1)
    return tuple(out)

def main():
    # Every signed basis product reduces to XOR; bilinearity proves the map.
    basis=[tuple(int(i==j) for i in range(16)) for j in range(16)]
    for a in basis:
        for b in basis:
            gate(tuple(t&1 for t in mul(a,b))==parity_product(a,b),'CD parity multiplicative basis map')
    # A weight-three parity unit at d4 is outside the extreme-Walsh family:
    # one-point affine functions have weights1,7,9,15, never3.
    outside=tuple(int(i in (0,1,2)) for i in range(16))
    gate(sum(outside)==3 and parity_product(outside,outside)==basis[0],'outside H4 odd unit')
    gate(dot(outside,outside)%2==sum(outside)%2==1,'odd norm augmentation')
    for m in range(1,65): test_shell(m,four(m))
    x,y=test_shell(7,(2,1,1,1))
    gate(not any(a*a+b*b+c*c==7 for a in range(3) for b in range(3) for c in range(3)), 'W shell14 excluded')
    gate(dot(x,x)==14 and any(y),'full S shell14 present')
    xp=tuple(t&1 for t in x)
    gate(sum(xp)%2==dot(x,x)%2==0,'singular integral element reduces to nonunit')
    gate(not any(parity_product(xp,xp)),'even augmentation square zero')
    print(json.dumps({'status':'PASS','gates':GATES,'m_range':[1,64], 'explicit_norm14_element':x,'explicit_norm14_annihilator':y,'scope':'integer constructions and embeddings; infinitude proved in SYNTHESIS.md'},sort_keys=True))
if __name__=='__main__': main()
