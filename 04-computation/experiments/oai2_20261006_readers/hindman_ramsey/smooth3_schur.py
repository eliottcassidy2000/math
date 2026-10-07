#!/usr/bin/env python3
"""Schur / Folkman / Hindman-type colourings for 3-smooth summands.
(1) primitive S-unit solutions x+y=z, S={2,3}, gcd=1  (Gersonides check, finite range)
(2) Liouville colouring of S has no monochromatic Schur triple inside S
(3) SAT: is there a 2-colouring of [1..N] with no monochromatic {a,b,a+b}, a,b 3-smooth
    (variants: a!=b ; a==b allowed); find least N where UNSAT if any.
"""
import sys, itertools, math
from pysat.solvers import Solver

def smooth3(N):
    out=[]; a=1
    while a<=N:
        b=a
        while b<=N:
            out.append(b); b*=3
        a*=2
    return sorted(out)

def Omega(n):
    c=0
    while n%2==0: n//=2; c+=1
    while n%3==0: n//=3; c+=1
    assert n==1
    return c

# (1) primitive solutions
LIM=10**30
S=smooth3(LIM); Sset=set(S)
prim=[]
for i,x in enumerate(S):
    for y in S[i:]:
        if x+y>LIM: break
        if math.gcd(x,y)==1 and (x+y) in Sset:
            prim.append((x,y,x+y))
print("(1) primitive x<=y, x+y=z all 3-smooth, gcd 1, z<=1e30:", prim)

# (2) Liouville check
bad=[(x,y,z) for (x,y,z) in prim if (Omega(x)%2)==(Omega(y)%2)==(Omega(z)%2)]
print("(2) Liouville-monochromatic primitive triples:", bad, "(scaling by g multiplies all signs by lambda(g))")
# direct check in range
S2=smooth3(10**12); S2set=set(S2)
cnt=0
for i,x in enumerate(S2):
    for y in S2[i:]:
        if x+y>10**12: break
        if (x+y) in S2set and (Omega(x)-Omega(y))%2==0 and (Omega(x)-Omega(x+y))%2==0:
            cnt+=1
print("    direct count of Liouville-mono Schur triples (x<=y) with x+y<=1e12 in S:", cnt)

def schur_sat(N, allow_equal, extra_dist=False):
    S=smooth3(N); Sset=set(S)
    trip=[]
    for i,a in enumerate(S):
        for b in (S[i:] if allow_equal else S[i+1:]):
            if a+b>N: break
            trip.append((a,b,a+b))
    verts=sorted({v for t in trip for v in t})
    idx={v:k+1 for k,v in enumerate(verts)}
    s=Solver(name='cadical195')
    for (a,b,c) in trip:
        A,B,C=idx[a],idx[b],idx[c]
        s.add_clause(list({A,B,C}))
        s.add_clause(list({-A,-B,-C}))
    sat=s.solve()
    model=s.get_model() if sat else None
    s.delete()
    return sat, len(trip), len(verts), (model, idx) if sat else None

for allow in (False, True):
    print(f"(3) allow_equal={allow}")
    for e in range(1,41):
        N=2**e
        sat,nt,nv,_=schur_sat(N,allow)
        if not sat:
            print(f"   N=2^{e}={N}: UNSAT  (triples {nt}, vertices {nv})")
            # find least N by bisection between 2^(e-1) and 2^e
            lo,hi=2**(e-1),N
            while hi-lo>1:
                mid=(lo+hi)//2
                if schur_sat(mid,allow)[0]: lo=mid
                else: hi=mid
            print(f"   least UNSAT N = {hi}")
            break
        if e%5==0: print(f"   N=2^{e}: SAT (triples {nt}, vertices {nv})")
    else:
        print("   SAT through N=2^40")
