"""Independent graph/gcd descent audit with multiplicative atlas and subset DP."""
from pathlib import Path
from math import gcd,prod
from fractions import Fraction as F
from itertools import combinations
import json,hashlib
gates=0
def check(ok,message):
    global gates
    gates+=1
    if not ok:raise RuntimeError(message)

def primes(n):
    marked=[False]*(n+1);out=[]
    for p in range(2,n+1):
        if not marked[p]:
            out.append(p)
            for m in range(p*p,n+1,p):marked[m]=True
    return out

def sums():
    bag={1}
    for p in primes(356):
        if p%3==2:
            bag|={old*p**j for old in tuple(bag) for j in (1,2) if old*p**j<=356}
    return bag-{1}

def rank(rows):
    A=[[F(x) for x in row] for row in rows];r=0
    for c in range(len(A[0])):
        k=next((i for i in range(r,len(A)) if A[i][c]),None)
        if k is None:continue
        A[k],A[r]=A[r],A[k];v=A[r][c];A[r]=[x/v for x in A[r]]
        for i in range(r+1,len(A)):
            v=A[i][c]
            if v:A[i]=[x-v*y for x,y in zip(A[i],A[r])]
        r+=1
    return r

def main():
    print('Independent connected-decoder gcd descent and sharp controls')
    path=Path(__file__).with_name('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json')
    raw=path.read_bytes()
    check(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','pinned incoming profiles')
    levels=json.loads(raw)['levels']
    profiles={int(d):{(c,tuple(gs)) for c,gs in item['profiles']} for d,item in levels.items()}
    S=sums();atlas={(a,k-a) for k in S for a in range(1,(k+1)//2) if gcd(a,k)==1}
    check(len(atlas)==5855,'multiplicative atlas count')
    coeff={x for pair in atlas for x in pair}
    check(coeff==set(range(1,356)),'all oriented coefficients')
    g7=set(levels['6']['gcds'])
    clocks={a*c for c in g7 for a in coeff}
    check(len(g7)==42 and max(g7)==90,'incoming seven clocks')
    check(len(clocks)==6121 and max(clocks)==31950,'six-clock universe')
    # Independent prime-valuation form of the elementary outgoing-edge identity.
    for alpha in range(9):
        for beta in range(alpha,11):
            for gamma in range(11):
                check(max(alpha-gamma,0)<=max(beta-gamma,0),'valuation domination')
    Q=91**6
    for k,U in ((4,(1,3,4,9)),(5,(1,3,4,9,10)),(6,(1,3,4,9,10,16))):
        D=90*355**(7-k)
        core=tuple(D*u for u in U)+tuple(90*355**j for j in range(7-k))+(2,4,6,12)
        K=max(core);g=2*Q*K+1;row=core+(g,3*g)
        check(len(set(row))==13 and sum(row)<=Q*Q,'positive distinct physical box')
        check(g>2*Q*K,'dominance excludes every bounded crossing')
        check(all(v%7 for v in row),'literal strict 1/7 witness')
        masks=[0]*(1<<13)
        for mask in range(1,1<<13):
            bit=mask & -mask
            masks[mask]=gcd(masks[mask-bit],row[bit.bit_length()-1])
        check(masks[-1]==1 and masks[(1<<11)-1]==2,'primitive full row and core scale2')
        observed={size:0 for size in range(7,13)}
        for mask,c in enumerate(masks):
            size=mask.bit_count()
            if size in observed:
                word=tuple(sorted(gcd(c,row[i]) for i in range(13) if not(mask>>i&1)))
                check((c,word) in profiles[13-size],'full inherited shadow signature')
                observed[size]=max(observed[size],c)
        check(tuple(observed.values())==(90,6,6,3,2,1),'all subset maxima')
        edges=[];relations=[]
        parent=list(range(13))
        def find(v):
            while parent[v]!=v:v=parent[v]
            return v
        for i,j in combinations(range(13),2):
            d=gcd(row[i],row[j])
            if (row[i]+row[j])//d in S:
                edges.append((i,j));parent[find(i)]=find(j)
                relation=[0]*13;relation[i]=row[j]//d;relation[j]=-row[i]//d
                relations.append(relation)
        groups={}
        for i in range(13):groups.setdefault(find(i),set()).add(i)
        check({frozenset(v) for v in groups.values()}=={frozenset(range(11)),frozenset((11,12))},'actual all-scale decoder components')
        check(rank(relations)==11,'literal rational decoder rank')
        for size in (4,5,6):
            values=[masks[mask] for mask in range(1<<11) if mask.bit_count()==size]
            check(max(values)<=90*355**(7-size),'every smaller core subset cap')
        check(masks[(1<<k)-1]==D,'sharp selected subset')
        # Direct exact divisibility along every boundary edge for every core subset.
        for mask in range(1,(1<<11)-1):
            boundary=[(i,j) if mask>>i&1 else (j,i) for i,j in edges
                      if j<11 and bool(mask>>i&1)!=bool(mask>>j&1)]
            check(bool(boundary),'every proper subset has outgoing actual edge')
            d=masks[mask]
            for i,j in boundary:
                loss=d//gcd(d,row[j]);a=row[i]//gcd(row[i],row[j])
                check(a%loss==0 and loss<=355,'exact boundary loss')
        print('k',k,'sharp gcd',D,'physical maximum',K,'subset maxima',tuple(observed.values()),'decoder rank11')
    print('six-clock count',len(clocks),'maximum',max(clocks))
    print('PASS',gates,'always-active gates')

if __name__=='__main__':main()
