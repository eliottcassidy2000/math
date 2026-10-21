#!/usr/bin/env python3
"""Exact independent controls for distinguished path and rooted channels.
No imports from the repository's counting implementations.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import factorial, gcd, comb
import argparse, json, hashlib
from pathlib import Path

CHECKS=0

def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)

@lru_cache(None)
def odd_partitions(n, minimum=1):
    if n == 0:
        return ((),)
    return tuple((v,)+tail for v in range(minimum,n+1,2)
                 for tail in odd_partitions(n-v,v))

def cycle_data(part):
    mult={v:part.count(v) for v in set(part)}
    E=sum((v-1)//2 for v in part)+sum(gcd(part[i],part[j]) for i in range(len(part)) for j in range(i))
    z=1
    for v,c in mult.items():
        z*=v**c*factorial(c)
    return E,len(part),z

@lru_cache(None)
def channel(n,a,b):
    out=sum((Fraction(2**(a*E+b*K),z) for E,K,z in map(cycle_data,odd_partitions(n))),Fraction())
    check(out.denominator==1, f'nonintegral F({n},{a},{b})')
    return out.numerator

def path_count(n,k):
    return 2**comb(k-1,2)*channel(n-k,1,k)

def converse_path_count(n,k):
    r,delta=divmod(k,2)
    s,epsilon=divmod(n-k,2)
    if delta and epsilon:
        return 0
    return 2**(r*(r-1+delta+epsilon))*channel(s,2,k+epsilon)

def perm_action(n,p):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    index={e:b for b,e in enumerate(pairs)}
    dest=[]
    flipped=0
    for i,j in pairs:
        a,b=p[i],p[j]
        target=index[tuple(sorted((a,b)))]
        dest.append(1<<target)
        if a>b:
            flipped ^= 1<<target
    # Two byte-table transforms make literal orbit census cheap.
    tables=[]
    for start in range(0,len(pairs),8):
        width=min(8,len(pairs)-start)
        table=[]
        for value in range(1<<width):
            v=0
            for j in range(width):
                if value>>j&1:
                    v |= dest[start+j]
            table.append(v)
        tables.append((start,(1<<width)-1,table))
    def apply(mask):
        out=flipped
        for start,partmask,table in tables:
            out ^= table[(mask>>start)&partmask]
        return out
    return apply

def constrained_orbits(n,k):
    """Literal labelled masks; quotient only unmarked vertices by S_(n-k)."""
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    pinned=sum(1<<b for b,e in enumerate(pairs) if e in [(i,i+1) for i in range(k-1)])
    full=(1<<len(pairs))-1
    unseen={mask for mask in range(full+1) if mask&pinned==pinned}
    actions=[perm_action(n,tuple(range(k))+p) for p in permutations(range(k,n))]
    reversal=perm_action(n,tuple(reversed(range(k)))+tuple(range(k,n)))
    ordinary=fixed=0
    while unseen:
        mask=min(unseen)
        orbit={f(mask) for f in actions}
        check(all(v&pinned==pinned for v in orbit),'path constraint not preserved')
        unseen.difference_update(orbit)
        ordinary+=1
        fixed += reversal(full^mask) in orbit
    return ordinary,fixed

def rooted_class_controls(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    full=(1<<len(pairs))-1
    ps=list(permutations(range(n)))
    acts=[perm_action(n,p) for p in ps]
    unseen=set(range(full+1))
    classes=selfconverse=0
    weighted=[Fraction() for _ in range(n+1)]
    while unseen:
        mask=min(unseen)
        images=[f(mask) for f in acts]
        orbit=set(images)
        unseen.difference_update(orbit)
        classes+=1
        autos=[p for p,v in zip(ps,images) if v==mask]
        antis=[p for p,v in zip(ps,images) if v==full^mask]
        # Independent DFS of every directed simple path in this class.
        adjacent=[[False]*n for _ in range(n)]
        for b,(i,j) in enumerate(pairs):
            adjacent[i][j]=bool(mask>>b&1)
            adjacent[j][i]=not adjacent[i][j]
        h=[0]*(n+1)
        def walk(last,used,length):
            h[length]+=1
            for v in range(n):
                if not used>>v&1 and adjacent[last][v]:
                    walk(v,used|(1<<v),length+1)
        for v in range(n):
            walk(v,1<<v,1)
        for k in range(1,n+1):
            weighted[k]+=Fraction(h[k],len(autos))
        if not antis:
            continue
        selfconverse+=1
        unseen_roots=set(range(n))
        root_orbits=[]
        while unseen_roots:
            v=min(unseen_roots)
            roots={p[v] for p in autos}
            unseen_roots-=roots
            root_orbits.append(roots)
        central=[roots for roots in root_orbits if any(any(p[v]==v for p in antis) for v in roots)]
        check(len(central)==n%2,f'central root-orbit mismatch n={n},mask={mask}')
        check(all(sum(p[v]==v for v in range(n))==n%2 for p in antis),'anti fixed-root count')
    for k in range(1,n+1):
        expected=Fraction(2**(comb(n,2)-k+1),factorial(n-k))
        check(weighted[k]==expected,f'weighted path mass n={n},k={k}')
        check((weighted[k]==path_count(n,k)) == (n-k<=2),'free-action threshold')
    return {'n':n,'ordinary':classes,'self_converse':selfconverse,
            'weighted_path_mass':{str(k):str(weighted[k]) for k in range(1,n+1)}}

def multiply_linear_roots(roots):
    # D(z)=product(1-root*z), coefficients in ascending order.
    out=[1]
    for root in roots:
        nxt=[0]*(len(out)+1)
        for i,a in enumerate(out):
            nxt[i]+=a
            nxt[i+1]-=root*a
        out=nxt
    return out

def recurrence_controls():
    table=[]
    for m in range(13):
        spectrum=sorted({len(part) for part in odd_partitions(m)})
        target=[0] if m==0 else list(range(1 if m%2 else 2,m+1,2))
        check(spectrum==target,f'spectrum missing m={m}')
        D=multiply_linear_roots([2**K for K in spectrum])
        seq=[channel(m,1,k) for k in range(1,2*len(D)+7)]
        for t in range(len(D)-1,len(seq)):
            check(sum(D[j]*seq[t-j] for j in range(len(D)))==0,f'ordinary recurrence m={m}')
        # Weights grouped by K give exact partial-fraction numerators.
        weights={K:Fraction() for K in spectrum}
        for part in odd_partitions(m):
            E,K,z=cycle_data(part)
            weights[K]+=Fraction(2**E,z)
        check(all(v>0 for v in weights.values()),'nonpositive partial-fraction weight')
        if m>=3:
            c=Fraction(m*(m-1)*(m-2),3)*Fraction(2)**(4-2*m)
            check(weights[m-2]/weights[m]==c,f'first symmetry correction m={m}')
        table.append({'m':m,'denominator':D,'weights':{str(K):str(v) for K,v in weights.items()},'normalized_prefix':seq[:5]})
        s,epsilon=divmod(m,2)
        spec2=sorted({len(part) for part in odd_partitions(s)})
        D2=multiply_linear_roots([4**K for K in spec2])
        for delta in (0,1):
            if delta and epsilon:
                continue
            r0=1-delta
            seq2=[channel(s,2,2*r+delta+epsilon) for r in range(r0,r0+2*len(D2)+7)]
            for t in range(len(D2)-1,len(seq2)):
                check(sum(D2[j]*seq2[t-j] for j in range(len(D2)))==0,f'anti recurrence m={m},delta={delta}')
    return table

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--receipt',default='receipt.json')
    args=parser.parse_args()
    rows=[]
    for n in range(1,7):
        for k in range(1,n+1):
            ordinary,fixed=constrained_orbits(n,k)
            check(ordinary==path_count(n,k),f'path mismatch {(n,k,ordinary,path_count(n,k))}')
            check(fixed==converse_path_count(n,k),f'converse mismatch {(n,k,fixed,converse_path_count(n,k))}')
            check((ordinary+fixed)%2==0,'merged count not integer')
            rows.append({'n':n,'k':k,'ordinary':ordinary,'converse_fixed':fixed,'merged':(ordinary+fixed)//2})
    rooted=[rooted_class_controls(n) for n in range(1,7)]
    # Hostiles distinguish pointwise root marking, reversal, weighting, and mixed cycle multiplicities.
    check(path_count(3,1)==4 and 3*channel(3,1,0)==6,'ordinary differentiation hostile')
    check(converse_path_count(3,1)==2 and converse_path_count(4,1)==0,'rooted anti parity')
    check(converse_path_count(3,2)==2,'path reversal not pointwise double rooting')
    E,K,z=cycle_data((3,5))
    correct_deficit=Fraction(2**E*8,z)
    erroneous_deficit=Fraction(2**E*8,3*5*factorial(2))
    check(correct_deficit==Fraction(128,15) and erroneous_deficit==Fraction(64,15),'THM-260 denominator hostile')
    rational=recurrence_controls()
    g3seq=[0]+[channel(3,1,k) for k in range(1,4)]
    g3den=[1,-10,16]
    g3num=[sum(g3den[j]*g3seq[i-j] for j in range(min(i,2)+1)) for i in range(4)]
    check(g3num==[0,12,-32,0],'G3 displayed rational numerator')
    receipt={'status':'ALL CHECKS PASSED','checks':CHECKS,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'literal_path_orbits':rows,'self_converse_class_controls':rooted,'normalized_rational_series':rational,
      'hostiles':{'rooted_order3':4,'n_times_unlabelled_order3':6,'pointwise_two_root_converse':0,'reversed_two_path_order3':2,
      'THM260_mixed_cycle_correct':str(correct_deficit),'THM260_displayed_general_formula':str(erroneous_deficit)}}
    Path(args.receipt).write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'checks':CHECKS,'literal_path_cells':len(rows),'rooted_control_classes':sum(r['ordinary'] for r in rooted),'receipt':args.receipt}))

if __name__=='__main__':
    main()
