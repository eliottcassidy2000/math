#!/usr/bin/env python3
"""Prime-power constructor. Compressed output is an ordered list of (word, repeat)."""
from argparse import ArgumentParser
from math import lcm,prod
from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path
_spec=spec_from_file_location('apex_core',Path(__file__).with_name('chanpak-all-prime-apex-check.py'))
_core=module_from_spec(_spec);_spec.loader.exec_module(_core)
TILES,O,U,Z,inv,projection,mul=(_core.TILES,_core.O,_core.U,_core.Z,_core.inv,_core.projection,_core.mul)

def identity():return (0,0,1,0)
def power(x,k,q):
    y=identity()
    while k:
        if k&1:y=mul(y,x,q)
        x=mul(x,x,q);k>>=1
    return y

def word_channel(w,lam,q,E):
    if not lam:return sum((1 if t>0 else -1) for t in w)%q
    o=(1,0,lam,0);i=(0,1,lam,0); Y=[]
    for w0 in TILES:
        z=identity()
        for t in w0:z=mul(z,o if t=="O" else i,q)
        Y.append(z)
    YY=[power(y,E-1,q) for y in Y]
    r=identity()
    for j in w:r=mul(r,Y[j-1] if j>0 else YY[-j-1],q)
    return r

def glast(w,lam,q,E):
    z=word_channel(w,lam,q,E)
    return z if lam==0 else z[-1]

def edge_words(p):
    lines=[Z[:]]
    for j in range(1,(p-1)//2):
        t=(1+pow(j,-1,p))%p
        a,b=projection(p,t),projection(p,pow(t,-1,p))
        lines.append(a+b+inv(a)+inv(b))
    return lines

def solve(A,B,q,p):
    M=[list(a)+[b] for a,b in zip(A,B)];d=len(B)
    for j in range(d):
        r=next(i for i in range(j,d) if M[i][j]%p)
        M[j],M[r]=M[r],M[j]
        M[j]=[x*pow(M[j][j],-1,q)%q for x in M[j]]
        for i in range(d):
            if i!=j:
                c=M[i][j]
                M[i]=[(a-c*b)%q for a,b in zip(M[i],M[j])]
    return [m[-1] for m in M]

def construct(targets,even):
    """targets maps (p,e) to vector. Output positive SLP plus common inverse exponent."""
    assert len({p for p,e in targets})==len(targets)
    E=lcm(*(p**e*(p-1)*p**(e-1) for p,e in targets))
    product=prod(p**e for p,e in targets)
    raw=[(O[:],1)] if even else []
    for (p,e),target in targets.items():
        q=p**e;d=(p-1)//2
        assert len(target)==d
        lines=edge_words(p);M=[];V=[]
        for s,want in enumerate(target,1):
            z=pow(s,-1,q)
            M.append([(glast(w,(1+z)%q,q,E)+glast(w,(1-z)%q,q,E))%q for w in lines])
            V.append((want*pow(s,-1-int(even),q)-1)*s*s%q)
        A=solve(M,V,q,p)
        isolate=product//q*pow(product//q,-1,q)%product
        raw.extend((w,k*isolate%product) for w,k in zip(lines,A) if k)
    if not raw:raw=[([1],E)]
    return raw,E

def checkout(raw,E,p,e):
    q=p**e
    n=sum(k*sum(len(TILES[abs(t)-1])*(1 if t>0 else E-1) for t in w) for w,k in raw)
    def channel(lam):
        if lam==0:return sum(k*sum(1 if j>0 else E-1 for j in w) for w,k in raw)%q
        z=identity()
        for w,k in raw:z=mul(z,power(word_channel(w,lam,q,E),k,q),q)
        return z[-1]
    ans=[]
    for s in range(1,(p+1)//2):
        a=pow(s,-1,q)
        ans.append(pow(s,n+1,q)*(1+a*a*(channel((1+a)%q)+channel((1-a)%q)))%q)
    return n+1,ans

def check():
    for p in (3,5,7,11,17,23):
        for e in (2,3):
            T=[(i*i+3*p**(e-1)-6)%p**e for i in range((p-1)//2)]
            for parity in (False,True):
                w,E=construct({(p,e):T},parity);n,got=checkout(w,E,p,e)
                assert got==T and (n%2==0)==parity
    print("24 prime-power constructor checks including 17 and 23, both parities")
    T={(3,3):[17],(5,2):[0,18],(7,3):[21,0,15],(11,2):[0,21,24,7,0]}
    for parity in (False,True):
        w,E=construct(T,parity)
        for (p,e),expected in T.items():
            n,got=checkout(w,E,p,e)
            assert got==expected and (n%2==0)==parity
    print("joint prime-power constructor (27,25,343,121) passed, both parities")

if __name__=="__main__":
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("--prime",type=int);ap.add_argument("--power",type=int,default=1)
    ap.add_argument("--even",action="store_true");ap.add_argument("target",type=int,nargs="*")
    a=ap.parse_args()
    if a.prime:
        raw,E=construct({(a.prime,a.power):a.target},a.even)
        for w,k in raw:
            print(k,[(TILES[abs(i)-1],1 if i>0 else E-1) for i in w])
    else:check()
