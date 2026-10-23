#!/usr/bin/env python3
"""Interval-chain product, strong connectivity, fixed scores, small genus2 collisions."""
from independent_audit import H,components,need
import itertools,math,json
from collections import Counter,defaultdict
from pathlib import Path
OUT=Path(__file__).resolve().parent

def compositions(n):
    if n==0:yield ()
    for p in range(2,n+1):
        for rest in compositions(n-p):yield (p,)+rest

def chain(parts):
    n=1+sum(parts);z=[[0 if i==j else (1 if i<j else -1) for j in range(n)] for i in range(n)];a=0
    for p in parts:b=a+p;z[a][b]=-1;z[b][a]=1;a=b
    return z

def c4(z):
    out=0
    for s in itertools.combinations(range(len(z)),4):
        for rest in itertools.permutations(s[1:]):
            word=(s[0],)+rest+(s[0],)
            out+=all(z[i][j]==1 for i,j in zip(word,word[1:]))
    return out

def rootpair(z):
    rr=list(itertools.permutations(range(len(z)),2));out=Counter()
    for a,b in rr:
        for c,d in rr:
            g=int(a==c)-int(a==d)-int(b==c)+int(b==d)
            w=z[a][c]-z[a][d]-z[b][c]+z[b][d]
            out[g,w]+=1
    return out

def main():
    tested=0;collisions=[];rows=[]
    for N in range(2,13):
        groups=defaultdict(list)
        for parts in compositions(N):
            z=chain(parts);actual=H(z);expected=math.prod(2**(p-1)+1 for p in parts)
            need(actual==expected,('interval product',parts,actual,expected))
            need(len(components(z))==1,('not strong',parts))
            scores=[sum(v==1 for v in row) for row in z]
            need(scores==[N-1]+[N-i for i in range(1,N)]+[1],('fixed labelled scores',parts,scores))
            groups[c4(z)].append((actual,parts));tested+=1
        for four,group in groups.items():
            if len({h for h,p in group})>1:
                a=group[0];b=next(x for x in group if x[0]!=a[0])
                need(rootpair(chain(a[1]))==rootpair(chain(b[1])),('genus2 candidate failed',a,b))
                collisions.append(dict(n=N+1,c4=four,H1=a[0],parts1=a[1],H2=b[0],parts2=b[1]))
        rows.append(dict(n=N+1,chains=sum(map(len,groups.values()))))
    primes=[p for p in range(5,100) if all(p%d for d in range(2,math.isqrt(p)+1))]
    b=[(2**p+1)//3 for p in primes]
    need(all(x>1 and x%3 for x in b),'bad prime quotient')
    need(all(math.gcd(x,y)==1 for x,y in itertools.combinations(b,2)),'prime quotients not coprime')
    replacements=[]
    for mask in range(4):
        parts=tuple(x for i,p in enumerate((5,7)) for x in ((p+1,) if mask>>i&1 else (2,)*((p+1)//2)))
        z=chain(parts);actual=H(z);need(actual==math.prod(2**(x-1)+1 for x in parts),'replacement product')
        replacements.append(dict(mask=mask,parts=parts,n=len(z),H=actual))
    need(len({x['H'] for x in replacements})==4,'replacement H collision')
    report=dict(status='PASS',chains_checked=tested,rows=rows,coprime_primes=primes,replacement_family=replacements,explicit_genus2_collisions=collisions)
    (OUT/'interval_chain_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
