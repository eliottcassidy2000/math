#!/usr/bin/env python3
import itertools as it, math, json, hashlib
from pathlib import Path
P=Path(__file__).resolve().parent
ROWS=((2,0),(0,2),(2,2),(2,-2),(1,2),(1,-2),(2,1),(2,-1))
OFF=((0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1))
def check(t,msg):
    if not t: raise RuntimeError(msg)
def C(n):
    g=lambda d:math.gcd(n,d)
    return 4*g(4)+4*g(6)+2*g(2)**2+2*g(5)+g(2)*g(4)-2-4*g(2)
def formula(n):return n*(n*n-(4*math.gcd(n,2)+4)*n+C(n))
def offsets(u,v,n):return tuple((a*u+b*v)%n for a,b in OFF)
def square_count(n):
    mask=sum(1<<q for q in {x*x%n for x in range(n)}); full=(1<<n)-1
    rotated=[((mask>>d)|(mask<<(n-d)))&full for d in range(n)]
    out=0
    for u,v in it.product(range(n),repeat=2):
        off=offsets(u,v,n)
        if len(set(off))!=9:continue
        valid=full
        for d in off:valid &=rotated[d]
        out+=valid.bit_count()
    return out

def run():
    agg={}; subset_rows=[]
    for k in range(2,9):
        for ids in it.combinations(range(8),k):
            rs=[ROWS[i] for i in ids]
            d1=math.gcd(*(x for r in rs for x in r))
            minors=[a*d-b*c for (a,b),(c,d) in it.combinations(rs,2)]
            d2=math.gcd(*minors)//d1
            check(d2>0,'rank two')
            key=(d1,d2);agg[key]=agg.get(key,0)+(-1)**k
            subset_rows.append([list(ids),d1,d2,(-1)**k])
    agg={k:v for k,v in agg.items() if v}
    expected={(1,4):4,(1,6):4,(2,2):2,(1,5):2,(2,4):1,(1,1):-2,(1,2):-4}
    check(agg==expected,'Smith aggregates')
    controls=[]
    for n in range(1,121):
        valid=0
        for u,v in it.product(range(n),repeat=2):
            distinct=len(set(offsets(u,v,n)))==9
            outside=all((a*u+b*v)%n for a,b in ROWS)
            check(distinct==outside,'collision union')
            valid+=distinct
        check(n*valid==formula(n),'ordinary exact formula')
        controls.append([n,n*valid])
    check(all(C(n)==C(n+60) for n in range(1,61)),'period60')
    witnesses={d:next(n for n in range(1,61) if C(n)!=C(n+d)) for d in (12,20,30)}
    # Independent bitset intersection of translated square sets, no recurrence used in census.
    square_controls=[]
    for p,K in ((2,9),(3,5)):
        for k in range(1,K+1):
            n=p**k;shift=3 if p==2 else 1
            predicted=sum(formula(p**(j-shift)) for j in range(k,0,-2) if j>=shift)
            actual=square_count(n)
            check(actual==predicted,'independent square-center census')
            square_controls.append([p,k,actual])
    gf=[]
    for p,shift,num in ((2,3,{7:1536}),(3,1,{3:216,4:5832})):
        coeff=[sum(formula(p**(j-shift)) for j in range(k,0,-2) if j>=shift) for k in range(51)]
        out=coeff[:]
        for d,c in ((2,1),(1,p),(1,p*p),(1,p**3)):
            out=[out[k]-(c*out[k-d] if k>=d else 0) for k in range(len(out))]
        check(all(x==num.get(k,0) for k,x in enumerate(out)),'GF numerator')
        gf.append({'p':p,'through':50,'first13':coeff[:13]})
    return {'status':'PASS','smith_aggregates':{str(k):v for k,v in sorted(agg.items())},'smith_subsets':subset_rows,'ordinary_controls':controls,'C_residues_0_to_59':[C(n) for n in range(60)],'proper_period_witnesses':witnesses,'square_controls':square_controls,'gf_controls':gf,'producer_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((P/'source-snapshots').iterdir())}}
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    out=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if a.write:(P/'certificate.json').write_text(out)
    else:check((P/'certificate.json').read_text()==out,'certificate mismatch')
    print('PASS: all-modulus Smith formula, collision equivalence, square-bitset census, generating functions, minimal period controls')
