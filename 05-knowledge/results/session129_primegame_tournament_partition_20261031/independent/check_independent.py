#!/usr/bin/env python3
"""Independent Session129 audit; no producer-code imports; stdlib only."""
import json, math, sys
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import permutations, combinations, product
from pathlib import Path

def require(c,s):
    if not c: raise RuntimeError(s)

@lru_cache(None)
def pairs(n): return tuple(combinations(range(n),2))
@lru_cache(None)
def orders(n): return tuple(permutations(range(n)))
def matrix(n, bits):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs(n)):
        a[i][j]=(bits>>k)&1; a[j][i]=1-a[i][j]
    return a

def code(a,p):
    return sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs(len(p))))

@lru_cache(None)
def orbit_map(n):
    result={}
    for b in range(1<<len(pairs(n))):
        if b in result: continue
        a=matrix(n,b)
        orbit={code(a,p) for p in orders(n)}
        root=min(orbit)
        result.update((x,root) for x in orbit)
    return result

def merged_key(n,b):
    m=orbit_map(n); flip=((1<<len(pairs(n)))-1)^b
    return min(m[b],m[flip])

def n4_h_key(a):
    return sum(all(a[p[i]][p[i+1]] for i in range(3)) for p in orders(4))

def full_labelled_audit():
    counts=defaultdict(Counter); pure=defaultdict(lambda:True)
    census={}
    for n in (3,4,5):
        roots=set(orbit_map(n).values())
        full=(1<<len(pairs(n)))-1
        sc=sum(orbit_map(n)[r^full]==r for r in roots)
        census[n]=(len(roots),sc,len({merged_key(n,r) for r in roots}))
    require(census=={3:(2,2,2),4:(4,2,3),5:(12,8,10)},'full-labelled census')
    n=5
    for b in range(1<<len(pairs(n))):
        a=matrix(n,b); root=merged_key(n,b)
        # H distinguishes the three merged n4 classes, avoiding producer canonical IDs.
        deletion_h={v:n4_h_key([[a[i][j] for j in range(n) if j!=v] for i in range(n) if i!=v]) for v in range(n)}
        for p in orders(n):
            if not all(a[p[i]][p[i+1]] for i in range(n-1)): continue
            counts[root][deletion_h[p[-1]]]+=1
            pure[root] &= all(a[p[i]][p[j]]==a[p[n-1-j]][p[n-1-i]] for i,j in pairs(n))
    rows=[]
    for root,c in sorted(counts.items()):
        require(all(c[h]%math.factorial(n)==0 for h in (1,3,5)),'label-normalization integrality')
        raw=tuple(c[h]//math.factorial(n) for h in (1,3,5))
        g=reduce(math.gcd,raw); v=tuple(x//g for x in raw)
        rows.append(dict(independent_root=root,face_vector=raw,primitive_face_vector=v,pure_blue=pure[root],label_7_17_55=7**v[0]*17**v[1]*55**v[2]))
    expected={(1,0,0),(1,1,0),(2,3,5),(0,1,0),(2,0,7),(1,4,4),(1,0,10),(0,2,3),(0,5,8),(0,0,1)}
    require({tuple(r['primitive_face_vector']) for r in rows}==expected,'n5 profile list')
    require({tuple(r['primitive_face_vector']) for r in rows if r['pure_blue']}=={(1,0,0),(0,1,0),(0,0,1)},'pure-blue axes')
    require({r['label_7_17_55'] for r in rows if r['pure_blue']}=={7,17,55},'chosen arithmetic labels')
    require(len({r['label_7_17_55'] for r in rows})==10,'arithmetic encoding injective')
    claimed=json.loads((Path(__file__).resolve().parents[1]/'tournaments/census.json').read_text())
    require({tuple(r['face_vector']) for r in rows}=={tuple(r['face_vector']) for r in claimed['5']['nodes']},'producer raw profiles disagree')
    # Order six: construct isomorphism orbits in the FULL 32768-labelled universe,
    # then count all Hamiltonian paths per ordinary representative / automorphism.
    # This is independent of the producer's 1024 fixed-path masks.
    n6map=orbit_map(6); n6roots=set(n6map.values()); full6=(1<<15)-1
    n6census=(len(n6roots),sum(n6map[r^full6]==r for r in n6roots),len({merged_key(6,r) for r in n6roots}))
    require(n6census==(56,12,34),'full labelled n6 census')
    n5order={tuple(r['primitive_face_vector']):r['node'] for r in claimed['5']['nodes']}
    coordinate_by_root={r['independent_root']:n5order[tuple(r['primitive_face_vector'])] for r in rows}
    six=defaultdict(lambda:[Fraction(0) for _ in range(10)])
    for r in sorted(n6roots):
        a=matrix(6,r); orbitcodes=[code(a,p) for p in orders(6)]
        aut=orbitcodes.count(r); upper=merged_key(6,r)
        lower={v:coordinate_by_root[merged_key(5,code(a,tuple(i for i in range(6) if i!=v)))] for v in range(6)}
        for p in orders(6):
            if all(a[p[i]][p[i+1]] for i in range(5)):
                six[upper][lower[p[-1]]]+=Fraction(1,aut)
    sixrows=[]
    for root,v in sorted(six.items()):
        require(all(x.denominator==1 for x in v),'n6 automorphism normalization integral')
        v=tuple(int(x) for x in v);g=reduce(math.gcd,v); prim=tuple(x//g for x in v)
        sixrows.append(dict(independent_root=root,face_vector=v,primitive_face_vector=prim))
    require(len({tuple(r['primitive_face_vector']) for r in sixrows})==34,'independent n6 primitive injectivity')
    require({tuple(r['face_vector']) for r in sixrows}=={tuple(r['face_vector']) for r in claimed['6']['nodes']},'producer n6 profiles disagree')
    return dict(census=census,n5_profiles=rows,independent_labelled_graphs=1024,
        hamiltonian_path_incidences=sum(sum(c.values()) for c in counts.values()),
        independent_n6_labelled_universe=32768,n6_census=n6census,n6_profiles=sixrows)

def rational_rank(rows):
    rows=[[Fraction(x) for x in row] for row in rows]
    rank=0
    for col in range(len(rows[0])):
        piv=next((i for i in range(rank,len(rows)) if rows[i][col]),None)
        if piv is None: continue
        rows[rank],rows[piv]=rows[piv],rows[rank]
        z=rows[rank][col]; rows[rank]=[x/z for x in rows[rank]]
        for i in range(len(rows)):
            if i!=rank and rows[i][col]:
                z=rows[i][col]; rows[i]=[x-z*y for x,y in zip(rows[i],rows[rank])]
        rank+=1
    return rank

def partition_audit(f):
    primes=(2,3,5,7,11,13,17,19,23,29)
    def val(n,p):
        e=0
        while n%p==0: n//=p;e+=1
        return e
    delta={k:tuple(val(v.numerator,p)-val(v.denominator,p) for p in primes) for k,v in zip('ABCDEFGHIJKLMN',f)}
    def inner(x,y):return sum(a*b for a,b in zip(x,y))
    phi=(-4,0,-4,2,2,2,2,3,-1,3)
    groups={'P':'JK','Q':'AIMN','R':'CEFH','S':'BDGL'}
    increments={k:{inner(phi,delta[c]) for c in g} for k,g in groups.items()}
    require(increments==dict(P={0},Q={-2},R={1},S={0}),'partition-constant character')
    constraints=[delta[c] for c in 'JKBDGL']
    for g in ('AIMN','CEFH'):
        constraints.extend(tuple(x-y for x,y in zip(delta[c],delta[g[0]])) for c in g[1:])
    require(rational_rank(constraints)==9,'character uniqueness rank')
    require(rational_rank(list(delta.values()))==10,'full increment rank')
    generators=('JK','ABDEFG','AABCEFH K L'.replace(' ',''),'ABEEFFLM')
    for w in generators:
        require(all(sum(delta[c][i] for c in w)==0 for i in range(10)),'kernel generator')
    # Symbolic completeness: all kernel vectors are uniquely determined by
    # nonnegative parameters (j,d,c,m); columns below give every coordinate.
    rows={k:tuple(w.count(k) for w in generators) for k in 'ABCDEFGHIJKLMN'}
    require(rows['J']==(1,0,0,0) and rows['D']==(0,1,0,0) and rows['C']==(0,0,1,0) and rows['M']==(0,0,0,1),'free nonnegative decoder')
    # Independent columns span rational kernel (four independent zero vectors).
    require(rational_rank(list(zip(*rows.values())))==4,'kernel generator rank')
    T=(1,0,1,0,0,0,0,0,1,0)
    require({c:inner(T,d) for c,d in delta.items()}=={c:int(c=='N') for c in delta},'N counter')
    def step(n):
        for c,g in zip('ABCDEFGHIJKLMN',f):
            if n%g.denominator==0:return c,n//g.denominator*g.numerator
    def trace(n,k):
        w=''
        for _ in range(k):c,n=step(n);w+=c
        return w,n
    require(trace(637,4)==('AIMN',55),'Q macro positive')
    require(trace(153,4)==('CHEF',539),'R macro positive')
    possible=[]
    for p in permutations('BDGL'):
        n=34
        for c in p:
            g=f[ord(c)-ord('A')]
            if n%g.denominator:break
            n=n//g.denominator*g.numerator
        else:possible.append(''.join(p))
    require(not possible,'BDGL all-permutation hostile')
    macro=[tuple(sum(delta[c][i] for c in groups[g]) for i in range(10)) for g in ('Q','R','S')]
    decoder=(T,tuple(int(p==11)-T[i] for i,p in enumerate(primes)),tuple(int(p==13)+T[i] for i,p in enumerate(primes)))
    require([[inner(d,m) for m in macro] for d in decoder]==[[1,0,0],[0,1,0],[0,0,1]],'integral left inverse')
    require(trace(4301,6)==('GBDEFA',4301),'genuine non-JK global cycle')
    samples=0
    for a in range(4):
      for u in range(1,4):
       for v in range(4):
        for w in range(4):
         for z in range(1,4):
          if a*w:continue
          n=2**a*11**u*13**v*17*19**w*23**z
          require(trace(n,6)==('GBDEFA',n),'GBDEFA exact guard family')
          samples+=1
    require(samples==252,'cycle family sample count')
    six_words=Counter();next_blocks=defaultdict(set)
    label={c:k for k,g in groups.items() for c in g}
    for x in product(range(3),repeat=10):
        n=math.prod(p**e for p,e in zip(primes,x))
        word,end=trace(n,6)
        next_blocks[label[word[0]]].add(label[word[1]])
        if end==n and set(word)==set('ABDEFG'):six_words[word]+=1
        guard=(x[1]==x[2]==x[3]==x[9]==0 and x[6]==1 and x[4]>=1 and x[8]>=1 and x[0]*x[7]==0)
        require((word=='GBDEFA')==guard,'all-capped exact six-cycle guard')
    require(set(six_words)=={'GBDEFA'[k:]+'GBDEFA'[:k] for k in range(6)},'saturated six-cycle word classification')
    require(dict(next_blocks)=={'P':{'P','Q','R'},'Q':{'P','Q','R','S'},'R':{'P','Q','R','S'},'S':{'P','Q','R','S'}},'global block graph saturation')
    return dict(full_increment_rank=10,character_constraint_rank=9,
                global_six_cycle='4301 --GBDEFA-->4301',cycle_family_samples=samples,block_codec_integral_left_inverse=True,
                saturated_six_cycle_words=dict(six_words),six_cycle_capped_universe=3**10,
                block_increments={k:sorted(v) for k,v in increments.items()},
                positive_priority_macros=['637 --AIMN-->55','153 --CHEF-->539'],
                product_integral_but_no_permutation_seed=34,kernel_generators=generators)

def encoding_audit(f):
    p=(2,3,5,7,11,13,17,19,23,29);q=(31,37,41,43,47,53,59,61,67,71)
    swap=dict(zip(p,q));swap.update(zip(q,p))
    def encode(n):
        result=1
        for a,b in swap.items():
            while n%a==0:n//=a;result*=b
        return result*n
    def step(n,rules):
        for i,g in enumerate(rules):
            if n%g.denominator==0:return i,n//g.denominator*g.numerator
    renamed=[Fraction(encode(g.numerator),encode(g.denominator)) for g in f]
    checks=0; supported=0
    for seed in range(1,501):
        residual=seed
        for prime in p:
            while residual%prime==0:residual//=prime
        supported+=residual==1
        n=seed;m=encode(n)
        for _ in range(100):
            i,n=step(n,f);j,m=step(m,renamed)
            require(i==j and m==encode(n),'whole-prime-permutation covariance')
            checks+=1
    require(supported==255 and checks==50000,'exact seed/step counts')
    require(encode(55)==1927 and math.isqrt(8*1927+1)**2!=8*1927+1,'triangular hostile')
    a,b=0,1
    while b<1927:a,b=b,a+b
    require(b!=1927,'Fibonacci hostile')
    require(step(2,f)==(11,15) and step(2,list(reversed(f)))==(0,110),'priority hostile')
    return dict(supported_seeds_to_500=supported,all_positive_seeds_to_500=500,
                permutation_conjugate_steps=checks,renaming_extension='swap P with Q; fix only primes outside P union Q')

def density_audit():
    p2=(3,7,13,17,19,23,29);p6=(3,5,7,17,29)
    base2=(11,13);base6=(4301,17765,16302,9867,8671,23023)
    rho2=math.prod(Fraction(p-1,p) for p in p2)
    rho6=math.prod(Fraction(p-1,p) for p in p6)*Fraction(37,38)
    d2=rho2*sum(Fraction(1,b) for b in base2)
    d6=rho6*sum(Fraction(1,b) for b in base6)
    require(d2==Fraction(2654208,36409529),'period2 density')
    require(d6==Fraction(37396469504,151504041213525),'period6 density')
    require(d2+d6==Fraction(144064009693952,1969552535775825),'total periodic-point density')
    points=set();counts={}
    for period,bases,ps in ((2,base2,p2),(6,base6,p6)):
        count=0
        for b in bases:
            for scale in range(1,100000//b+1):
                if math.gcd(scale,math.prod(ps))!=1 or (period==6 and scale%38==0):continue
                require(b*scale not in points,'phase-family disjointness')
                points.add(b*scale);count+=1
        counts[period]=count
    require(counts=={2:7288,6:28},'direct scale parameter count')
    return dict(period2=str(d2),period6=str(d6),total=str(d2+d6),direct_counts_to_100000=counts)

def basic_arithmetic():
    f=[Fraction(a,b) for a,b in [(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]]
    names=dict(zip('ABCDEFGHIJKLMN',f))
    blocks={k:math.prod(names[x] for x in word) for k,word in [('P','JK'),('Q','AIMN'),('R','CEFH'),('S','BDGL')]}
    require(blocks==dict(P=Fraction(1),Q=Fraction(55,637),R=Fraction(5929,1683),S=Fraction(585,34)),'four literal block products')
    # Harmonic reciprocal identity and generalized-Fibonacci/Cassini boundary.
    for n in range(1,1001):
        w=(1,n,n+1,2*n+1)
        require(w[0]+w[1]==w[2] and w[1]+w[2]==w[3],'recurrence window')
        a,b,c=n,n+1,n*(n+1)
        require(Fraction(1,a)==Fraction(1,b)+Fraction(1,c),'harmonic identity')
        require((math.gcd(a,b),math.gcd(a,c),math.gcd(b,c))==(1,n,n+1),'gcd decoder')
        require((abs(1+n-n*n)==1)==(n in (1,2)),'unit Cassini boundary')
    return dict(block_products={k:str(v) for k,v in blocks.items()},harmonic_states_checked=1000,partition=partition_audit(f),encoding=encoding_audit(f),density=density_audit())

if __name__=='__main__':
    if '--test-failing-gate' in sys.argv: require(False,'intentional active-gate negative control')
    print(json.dumps(dict(status='passed',tournaments=full_labelled_audit(),arithmetic=basic_arithmetic()),indent=2,sort_keys=True))
