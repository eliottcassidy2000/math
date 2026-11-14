"""Target-first construction for a complete line at the maximum outer depth."""
from functools import cache
from math import gcd
from physical_seed import decode as terrace_decode, adjacency as terrace_rows
from physical_prefix import row
from bridge import bridge, h_f2, rooted_profile

# Arbitrary core in the first row is not a terrace replacement.
CORE=(2,24,3,53,37,7)
DEFAULT=(CORE,(2,4,1))

def ordsum(parts):
    P=[tuple(r) for r in parts]
    n=sum(map(len,P))
    out=[]
    pos=0
    for R in P:
        later=((1<<n)-1)^((1<<(pos+len(R)))-1)
        out += [(r<<pos)|later for r in R]
        pos+=len(R)
    return out

def actual(t):
    n=len(t)
    return all(not(t[i]>>i&1) for i in range(n)) and all(
        (t[i]>>j&1)^(t[j]>>i&1) for i in range(n) for j in range(i+1,n))

def strong(t):
    n=len(t)
    for q in range(n):
        seen=1<<q
        U=seen
        while U:
            c=U&-U;U^=c
            nb=t[c.bit_length()-1]&~seen
            seen|=nb
            U|=nb
        if seen!=(1<<n)-1:return False
    return True

def pair_compiler(bank):
    bank=tuple(tuple(x) for x in bank)
    assert all(actual(R) and strong(R) and len(R)>=3 for R in bank)
    menus=tuple((R,*h_f2(R)) for R in bank)

    @cache
    def reservoirs(a,L):
        if a==1:
            if L>=2 and not L&(L-1):
                return (),L.bit_length()-1
            return None
        for core,h,f in menus:
            if a%h==0 and L%f==0:
                z=reservoirs(a//h,L//f)
                if z is not None:
                    C,r=z
                    return (core,)+C,r
        return None

    @cache
    def marked_R(R,v):
        return rooted_profile(R,v)[1]

    def rooted(a,b):
        # At a prime cofactor gcd(a,b)=1, so there is no nontrivial
        # strong component before the marked component.
        if gcd(a,b)!=1:return None
        # Route through a singleton root, then try marking an actual strong core.
        for B,h in (((0,),1),)+tuple((R,h) for R,h,_ in menus if a%h==0):
            for v in range(len(B)):
                Q=marked_R(B,v)
                if b%Q:continue
                z=reservoirs(a//h,b//Q)
                if z is None:continue
                C,r=z
                T=ordsum((B,)+C)
                P=bridge(T,v,r)
                return P,{"marked_component":B,"root":v,"forward_components":C,"bridge_width":r}
        return None
    return rooted

def target_line(p):
    if p<159*159 or not p&1:return None
    k=0
    u,v=row(0)
    while True:
        U,V=u+4*v,2*u+4*v
        if (U+2*V)**2>p:break
        u,v=U,V;k+=1
    a0=p*pow(u,-1,2*v)%(2*v)
    amax=p//(u+2*v)
    pairs=[(a,(p-u*a)//v) for a in range(a0,amax+1,2*v)]
    assert len(pairs)<=49
    return k,pairs

def enclosed_prefix(P,k):
    P=list(P)
    for q,r in ((1,1),)*k+((7,2),):
        for i in range(q):P=bridge(P,len(P)-1,0)
        P=bridge(P,len(P)-1,r)
    return P

def construct(p,bank=DEFAULT,prefer_bridge=False):
    z=target_line(p)
    if not z:return None
    k,points=z
    builders=pair_compiler(bank)
    for a,b in points:
        y=None
        if prefer_bridge:y=builders(a,b)
        if y is None:
            w=terrace_decode(a,b)
            if w is not None:y=(terrace_rows(w),{"gap_word":w})
        if y is None and not prefer_bridge:y=builders(a,b)
        if y is None:continue
        P,proof=y
        final=enclosed_prefix(P,k)
        assert actual(final) and strong(final)
        return {"cofactor":p,"H":7*p,"depth":k,"state":(a,b),
                "source":proof,"prefix_rows":P,"witness_rows":final,
                "point_index":points.index((a,b))}
    return None

if __name__=='__main__':
    from argparse import ArgumentParser
    from json import dumps
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("cofactor",type=int,nargs="*", default=(1137067,1512281))
    for p in ap.parse_args().cofactor:
        D=construct(p)
        print(dumps(D))
        if p in (1137067,1512281):
            state=(1105,2634) if p==1137067 else (23,4736)
            assert tuple(D["state"])==state
            assert rooted_profile(D["prefix_rows"],len(D["prefix_rows"])-1)==(
                sum(state),sum(state),-state[0])
