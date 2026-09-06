"""Independent literal half-lift audit of the ranked component theorem."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd, isqrt
import sys
sys.stdout.reconfigure(newline='\n')
N=0

def need(x, why):
    global N
    N+=1
    if not x: raise ArithmeticError(why)

def norm(x):
    x%=1
    return min(x,1-x)

def bad(y,p,q):
    return all(any(norm(w*(y+j)/2)<F(1,14) for w in (p,q)) for j in (0,1))

def raw_lengths(p,q):
    cuts={F(0),F(1)}
    for w in (p,q):
        for k in range(w+1):
            for e in (-1,1):
                z=F(7*k+e,7*w)
                if 0<z<1: cuts.add(z)
    cuts=sorted(cuts)
    arcs=[]
    for left,right in zip(cuts,cuts[1:]):
        if bad((left+right)/2,p,q):
            if arcs and arcs[-1][1]==left and bad(left,p,q):
                arcs[-1]=(arcs[-1][0],right)
            else: arcs.append((left,right))
    need(not bad(F(0),p,q) and not bad(F(1),p,q),'no circular endpoint merge')
    for l,r in arcs:
        need(l<r and not bad(l,p,q) and not bad(r,p,q),'strict component ends')
    return tuple(sorted((r-l for l,r in arcs),reverse=True))

def partial(v,k): return sum(v[:k],F(0))
def leq(u,v):
    return all(partial(u,k)<=partial(v,k) for k in range(1,max(len(u),len(v))+1))
def escape(u,v):
    return any(0<partial(u,k)>=partial(v,k) for k in range(1,len(u)+1))

FRONT={
 (7,13):(F(1,49),)*2,
 (19,25):(F(37,3325),)*2+(F(23,3325),)*2+(F(9,3325),)*2,
 (5,41):(F(2,287),)*6,
 (5,53):(F(2,371),)*6+(F(9,1855),)*2,
 (1,67):(F(2,469),)*10,
}

def main():
    sums={1}
    for p in range(2,135):
        if p%3==2 and all(p%d for d in range(2,isqrt(p)+1)):
            sums={s*p**e for s in sums for e in range(3) if s*p**e<=134}
    pairs=[(p,q) for q in range(1,68,2) for p in range(1,q,2)
           if p%3 and q%3 and gcd(p,q)==1 and p+q in sums]
    need(len(pairs)==46,'independent finite-head filters')
    bank={pair:raw_lengths(*pair) for pair in pairs}
    maxima={pair for pair,u in bank.items() if not any(u!=v and leq(u,v) for v in bank.values())}
    need(maxima==set(FRONT),'full weak-majorization frontier')
    for pair,v in FRONT.items(): need(bank[pair]==v,'complete profile')
    for pair,u in bank.items():
        need(any(leq(u,v) for v in FRONT.values()),'head covered')
        for g in (3,5):
            scaled=raw_lengths(g*pair[0],g*pair[1])
            expected=tuple(x/g for x in u for _ in range(g))
            need(scaled==expected,'literal physical odd-scale profile')
            need(leq(scaled,u),'splitting is weakly majorized')
    print('HEAD',len(pairs),'complete literal half-lift profiles; g=1,3,5')
    for pair,u in FRONT.items(): print('FRONT',pair,tuple(map(str,u)))

    # A different exact alphabet and lengths than the producer's paired bank.
    values=[F(0),F(1,1000),F(1,100),F(9,1855),F(2,371),F(23,3325),F(37,3325)]
    cases=0
    for half in combinations_with_replacement(values,4):
        u=tuple(sorted((x for x in half for _ in (0,1) if x),reverse=True))
        for pair,v in FRONT.items():
            reduced=sum(u)>=sum(v) or max(u,default=0)>=v[0] or (pair==(19,25) and partial(u,4)>=F(120,3325))
            need(escape(u,v)==reduced,'paired breakpoint iff')
        cases+=1
    print('PAIRED',cases,'fresh exact lists, eight entries with zero padding')

    u=(F(1,100),)*4+(F(1,2000),)*2
    need(not (sum(u)>=sum(FRONT[19,25]) or max(u)>=FRONT[19,25][0]),'abstract old miss')
    need(all(escape(u,v) for v in FRONT.values()),'abstract all-profile gain')
    need(partial(u,4)==F(1,25)>F(120,3325),'strict intermediate gain')
    need(leq((F(2,5),F(2,5)),(F(9,10),)),'many-to-one containment permitted')
    need(escape((F(9,10),),(F(9,10),)),'inclusive equality escape')
    need(not escape((),(F(9,10),)),'isolated-only geometry not certified')
    print('SCOPE: proof Sections1-4 accepted; abstract gain only; large body scout not replayed')
    print('PASS',N,'independent always-active gates')

if __name__=='__main__':main()
