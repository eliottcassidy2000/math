"""Exact Euclid-offset character syndromes and orientation controls."""
from math import isqrt,prod
from itertools import combinations
from pathlib import Path
import hashlib,json

def require(x,msg):
    if not x:raise ValueError(msg)

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def leg(a,p):
    require(prime(p) and p%2,'odd prime modulus')
    r=pow(a%p,(p-1)//2,p)
    return -1 if r==p-1 else r

def rank2(rows):
    r=0;rows=[list(x) for x in rows]
    for c in range(len(rows[0]) if rows else 0):
        pivot=next((i for i in range(r,len(rows)) if rows[i][c]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        for i in range(len(rows)):
            if i!=r and rows[i][c]:rows[i]=[a^b for a,b in zip(rows[i],rows[r])]
        r+=1
    return r

def check():
    ss=[2]
    for _ in range(6):ss.append(ss[-1]**2-ss[-1]+1)
    factors=[{2:1},{3:1},{7:1},{43:1},{13:1,139:1},{3263443:1},{547:1,607:1,1033:1,31051:1}]
    syl=[]
    for i,s in enumerate(ss):
        require(prod(q**e for q,e in factors[i].items())==s,'factor product')
        require(all(prime(q) for q in factors[i]),'prime factors')
        if i:
            require(s==prod(ss[:i])+1,'Euclid offset one')
            old=[p for d in factors[:i] for p in d if p%2]
            for p in old:
                require(s%p==1 and prod(leg(q,p)**e for q,e in factors[i].items())==1,'zero syndrome')
        if i>=2:require(s%4==3,'Sylvester residue three')
        syl.append({'index':i,'value':s,'factors':factors[i]})
    ps=[x for x in ss if x%4==3 and prime(x)]
    require(ps==[3,7,43,3263443],'prime-term scope')
    for p,q in combinations(ps,2):require(leg(p,q)==-1 and leg(q,p)==1,'transitive term orientation')
    cycles=sum(leg(a,b)==leg(b,c)==leg(c,a) for a,b,c in combinations(ps,3))
    require(cycles==0,'no cyclic triples')
    incoming=[[leg(q,p) for q in (13,139)] for p in (3,7,43)]
    outgoing=[[leg(p,q) for q in (13,139)] for p in (3,7,43)]
    require(all(a==b for a,b in incoming),'same incoming columns')
    require(all(a==-b for a,b in outgoing),'opposite root-character columns')
    require(rank2([[int(x==-1) for x in row] for row in incoming])==1,'rank one factor signature')
    ff=[2**(2**i)+1 for i in range(6)]
    ffactors=[{x:1} for x in ff[:5]]+[{641:1,6700417:1}]
    ferm=[]
    for j,f in enumerate(ff):
        require(prod(q**e for q,e in ffactors[j].items())==f,'Fermat factor product')
        require(all(prime(q) for q in ffactors[j]),'Fermat factor primality')
        if j:
            require(f==prod(ff[:j])+2,'Euclid offset two')
            old=[p for d in ffactors[:j] for p in d]
            for p in old:
                require(f%p==2 and prod(leg(q,p)**e for q,e in ffactors[j].items())==leg(2,p),'Fermat offset syndrome')
        ferm.append({'index':j,'value':f,'factors':ffactors[j]})
    for p,q in combinations(ff[2:5],2):require(leg(p,q)==1 and leg(q,p)==1,'Fermat graph symmetric clique not tournament')
    return {'status':'PASS','sylvester':syl,'fermat':ferm,'prime_term_tournament':ps,'cyclic_triples':cycles,
            'roots':[3,7,43],'factor_ports':[13,139],'incoming_factor_signs':incoming,'outgoing_root_signs':outgoing,
            'fermat_f5_syndrome':[leg(2,p) for p in ff[:5]],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Exact displayed factors and sign identities; all-order statements proved separately; no infinite Sylvester-prime claim.'}
if __name__=='__main__':print(json.dumps(check(),indent=2,sort_keys=True))
