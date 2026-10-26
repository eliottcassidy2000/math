"""Per-n analysis of supplied finite class tables: C3-invariance tests of H/|Aut|, labelled C3-structure counts, fibre spectrum,
Hartley information E[log2 N], labelled prime probability, THM-1430 tiling identity, tree-shape table."""
import os
DATA=os.environ.get("DATA",".")
import csv, math, sys, collections
from fractions import Fraction
for n in range(3, 10):
    rows=list(csv.DictReader(open(os.path.join(DATA,f't{n}.tsv')),delimiter='\t'))
    nf=math.factorial(n); m=n*(n-1)//2
    # THM-1430: sum H/|Aut| = 2^{C(n-1,2)}
    tiling=sum(Fraction(int(r['H']),int(r['aut'])) for r in rows)
    ok1430 = (tiling==2**((n-1)*(n-2)//2))
    # Finite tests only; |Aut| first fails at n=12. SCC field records sizes, not types.
    groups=collections.defaultdict(list)
    for r in rows: groups[r['c3hash']].append(r)
    bad=0
    for h,g in groups.items():
        keys=set((r['H'],r['aut'],r['c3'],tuple(sorted(map(int,r['compsizes'].split(','))))) for r in g)
        if len(keys)>1: bad+=1
    # class sizes per C3 iso-class: number of tournament iso classes sharing the hypergraph
    gsizes=collections.Counter(len(g) for g in groups.values())
    # labelled counts
    tot=sum(Fraction(nf,int(r['aut'])) for r in rows); assert tot==2**m
    labc3=sum(Fraction(nf,int(r['aut']))/round(2**float(r['logN'])) for r in rows)
    elog=sum(Fraction(nf,int(r['aut']))*Fraction(float(r['logN'])).limit_denominator(10**9) for r in rows)/2**m
    pprime=sum(Fraction(nf,int(r['aut'])) for r in rows if r['ncomp']=='1' and r['prime']=='1' and r['lin']=='0' and r['sig'].count('.')==n and r['sig'].startswith(f'P{n}:'))/2**m
    fib=collections.Counter(round(2**float(r['logN'])) for r in rows)
    shapes=collections.Counter((int(r['lin']),int(r['prime'])) for r in rows)
    nprime=sum(1 for r in rows if r['sig'].startswith(f'P{n}:') and r['sig'].count('.')==n)
    nprime_sc=sum(1 for r in rows if r['sig'].startswith(f'P{n}:') and r['sig'].count('.')==n and r['sc']=='1')
    print(f'n={n}: classes={len(rows)} THM1430={ok1430} C3-invariance violations(H,Aut,c3,SCC multiset)={bad} C3-iso classes={len(groups)}')
    print(f'   tournament classes per C3-class: {dict(sorted(gsizes.items()))}')
    print(f'   labelled C3-structures={labc3} ({float(labc3)/2**m:.4f} of 2^m) ; E[log2 N]={float(elog):.4f} ; P(prime, labelled)={float(pprime):.4f}; prime classes={nprime} (self-converse {nprime_sc})')
    print(f'   fibre spectrum N(T)->#classes: {dict(sorted(fib.items()))}')
    print(f'   tree shapes (#linear,#prime)->#classes: {dict(sorted(shapes.items()))}')
