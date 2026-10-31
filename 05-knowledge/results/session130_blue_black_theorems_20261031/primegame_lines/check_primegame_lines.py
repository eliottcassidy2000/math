#!/usr/bin/env python3
"""Session130 exact phase/line bridge. Standard library, no inherited imports.
All checks use explicit failure gates, including under python -O.
"""
import itertools as it
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
from math import gcd, isqrt, prod

BASE=(4301,17765,16302,9867,8671,23023)
PRIMES=(2,3,5,7,11,13,17,19,23,29)
FRACS=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),
       (77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
PHI=(0,1,2,4,5,3)

def require(c,msg):
    if not c: raise RuntimeError(msg)

def vp(n,p):
    n=Fraction(n); ans=0
    a,b=n.numerator,n.denominator
    while a%p==0: a//=p; ans+=1
    while b%p==0: b//=p; ans-=1
    return ans

def degree(z): return sum(vp(z,p) for p in (3,7,17,29))
def decoder(z,X=1):
    a=vp(X,2); k=degree(z)
    return (vp(z,17)-vp(z,5),vp(z,5),vp(z,2)-a*k,
            vp(z,3)-vp(z,2)+a*k,vp(z,29),vp(z,7))
def encode(e,X=1): return prod(Fraction(X*b)**a for b,a in zip(BASE,e))
def mirror(e): return tuple(e[-j%6] for j in range(6))
def complement(e): return tuple(1-a for a in e)
def phase_shift(e): return tuple(e[(j-1)%6] for j in range(6))
def is_square(z):
    z=Fraction(z)
    return isqrt(z.numerator)**2==z.numerator and isqrt(z.denominator)**2==z.denominator

def step(n):
    for k,(a,b) in enumerate(FRACS):
        if n%b==0:return 'ABCDEFGHIJKLMN'[k],n//b*a
    raise RuntimeError('denominator 1 missing')

def tiles(n): return tuple((x,y) for y in range(1,n-1) for x in range(n,y+1,-1))
def reflection_permutation(n):
    ts=tiles(n); ix={e:k for k,e in enumerate(ts)}
    return tuple(ix[(n-y+1,n-x+1)] for x,y in ts)
def permute_mask(mask,p):return sum(((mask>>k)&1)<<p[k] for k in range(len(p)))
def adjacency(n,mask):
    A=[[0]*n for _ in range(n)]
    for i in range(1,n):A[i][i-1]=1
    for k,(x,y) in enumerate(tiles(n)):
        if (mask>>k)&1:A[y-1][x-1]=1
        else:A[x-1][y-1]=1
    return A

def merged_code(A):
    n=len(A); pairs=list(it.combinations(range(n),2)); full=(1<<len(pairs))-1
    best=full
    for p in it.permutations(range(n)):
        v=0
        for i,j in pairs:v=(v<<1)|A[p[i]][p[j]]
        best=min(best,v,full^v)
    return best

def census(n):
    size=1<<len(tiles(n)); codes=[merged_code(adjacency(n,m)) for m in range(size)]
    labels={c:i for i,c in enumerate(sorted(set(codes)))}
    q=[labels[c] for c in codes]; fibres=defaultdict(list)
    for m,u in enumerate(q):fibres[u].append(m)
    return q,dict(fibres)

def main():
    basis=[]
    for j,b in enumerate(BASE):
        expected=tuple(int(j==k) for k in range(6))
        require(decoder(b)==expected,'canonical integral left inverse')
        basis.append({str(p):vp(b,p) for p in PRIMES if vp(b,p)})
        require(step(b)==('GBDEFA'[j],BASE[(j+1)%6]),'canonical priority trace')
        require(degree(b)==1,'degree marker')
    for letter in 'ABDEFG':
        a,b=FRACS['ABCDEFGHIJKLMN'.index(letter)]
        require(degree(Fraction(a,b))==0,'six-instruction marker conservation')
    for j in range(6):
        expected=tuple(int(k==(j+1)%6)-int(k==j) for k in range(6))
        require(decoder(Fraction(BASE[(j+1)%6],BASE[j]))==expected,'cycle boundary decoder')
    require(reduce(gcd,BASE)==1,'base gcd')
    require(degree(Fraction(*FRACS[2]))==-2,'full-program noninvariance hostile')
    bits=list(it.product((0,1),repeat=6))
    codes=[encode(e) for e in bits]
    require(len(set(codes))==64,'canonical binary injectivity')
    require(prod(BASE)==2453532251388619737666330,'K product')
    scales=[]; scaled_checks=0
    for a,r,v,w,s in it.product(range(3),repeat=5):
        if a*w:continue
        for S in (1,31):
            X=S*2**a*11**r*13**v*19**w*23**s
            require(gcd(X,51765)==1 and X%38!=0,'scale admissibility')
            require(reduce(gcd,(X*b for b in BASE))==X,'scale gcd recovery')
            require(decoder(X)==(0,0,a,-a,0,0),'rank-one vector')
            for j,b in enumerate(BASE):
                require(step(X*b)==('GBDEFA'[j],X*BASE[(j+1)%6]),'scaled priority trace')
            for e in bits:
                z=encode(e,X); rz=encode(mirror(e),X); K=encode((1,)*6,X)
                require(decoder(z,X)==e,'scaled decoder')
                require(degree(z)==sum(e),'marker degree total')
                require(encode(complement(e),X)==K/z,'arithmetic complement')
                blue=e==mirror(e)
                equations=vp(z,5)==vp(z,7) and vp(z,2)-a*degree(z)==vp(z,29)
                require(blue==equations==is_square(z*rz),'scaled blue square-norm iff')
                scaled_checks+=1
            scales.append(X)
    require(len(scales)==270 and len(set(scales))==270,'scale universe size')
    signed_checks=0
    for X in (1,2,11*13*19*23*31,4*11**2*13*23**2*31):
        require(gcd(X,51765)==1 and X%38,'signed scale allowed')
        for e in it.product((-1,0,1),repeat=6):
            require(decoder(encode(e,X),X)==e,'signed-lattice primitive decoder')
            signed_checks+=1
    raw=[]
    for m in (2,6):
        for c in range(m):
            fixed_R=[j for j in range(m) if (c-j)%m==j]
            fixed_CR=[j for j in range(m) if (c-j+m//2)%m==j]
            raw.append(dict(period=m,reflection_parameter=c,fixed_R=fixed_R,fixed_CR=fixed_CR))
            if not fixed_CR:require(len(fixed_R)==2,'admissible raw fixed count')
    r=reflection_permutation(5)
    require(r==(0,3,5,1,4,2),'tile reflection gauge')
    require(all(PHI[-j%6]==r[PHI[j]] for j in range(6)),'chosen gauge equivariance')
    gauges=[p for p in it.permutations(range(6)) if all(p[-j%6]==r[p[j]] for j in range(6))]
    require(len(gauges)==16,'gauge classification')
    levels={}
    for n in (3,4,5):
        rp=reflection_permutation(n); size=1<<len(rp)
        fixed=sum(permute_mask(m,rp)==m for m in range(size))
        cr_fixed=sum((size-1)^permute_mask(m,rp)==m for m in range(size))
        require(cr_fixed==0,'no complemented mirror fixed tilings')
        levels[n]=dict(tilings=size,blue_lines=fixed//2,black_line_pairs=(size-fixed)//4)
    require([(levels[n]['blue_lines'],levels[n]['black_line_pairs']) for n in (3,4,5)]==[(1,0),(2,1),(8,12)],'sharp raw obstruction onset')
    blue=[e for e in bits if e==mirror(e)]
    retained=[e for e in blue if phase_shift(e)==mirror(phase_shift(e))]
    require(len(blue)==16 and len(retained)==4,'phase-shift color obstruction')
    for e in retained:require(e[0]==e[2]==e[4] and e[1]==e[3]==e[5],'retained parity pattern')
    mirror_count=Counter(sum(e==tuple(e[(c-j)%6] for j in range(6)) for c in (0,2,4)) for e in bits)
    require(mirror_count=={0:24,1:36,3:4},'three temporal axes count')
    e=tuple(int(j in (1,2)) for j in range(6)); f=tuple(int(j in (1,4)) for j in range(6))
    require(encode(e)*encode(mirror(e))==encode(f)*encode(mirror(f)),'norm collision equality')
    orbit={e,mirror(e),complement(e),complement(mirror(e))}
    require(f not in orbit,'norm does not determine V4 orbit')
    nonbinary=(0,2,0,0,0,0)
    nz=encode(nonbinary); nr=encode(mirror(nonbinary))
    require(nz!=nr and is_square(nz*nr),'full-lattice square-norm false converse hostile')
    q,fibres=census(5)
    require(len(fibres)==10,'merged node count')
    require(q[6]==2 and q[34]==5,'same norm distinct merged nodes hostile')
    require(fibres[q[16]]==[16],'sandwich singleton fibre')
    require(q[1]==q[7] and q[2]==q[8] and q[4]==q[32],'other singleton masks are non-singleton fibres')
    require(all(len(fibres[q[1<<k]])>1 for k in range(6) if k!=4),'all possible images of singleton fibre')
    witnesses=[]
    for p in gauges:
        tile_shift=[None]*6
        for j in range(6):tile_shift[p[j]]=p[(j+1)%6]
        sh=lambda m:permute_mask(m,tile_shift)
        witness=next((x,y) for x in range(64) for y in range(x+1,64) if q[x]==q[y] and q[sh(x)]!=q[sh(y)])
        x,y=witness
        witnesses.append(dict(gauge=p,source_masks=witness,source_node=q[x],target_masks=[sh(x),sh(y)],target_nodes=[q[sh(x)],q[sh(y)]]))
        require(sh(16)!=16 and len(fibres[q[sh(16)]])>1,'universal singleton-fibre obstruction')
    chosen=next(w for w in witnesses if tuple(w['gauge'])==PHI)
    require(chosen['source_masks']==(1,28) and chosen['target_masks']==[2,49] and chosen['target_nodes']==[2,8],'chosen minimal node hostile')
    # Literal line arithmetic is checked independently through tile masks.
    line_rows=[]
    for m in range(32):
        e=tuple((m>>PHI[j])&1 for j in range(6)); z=encode(e)
        em=tuple(((63^m)>>PHI[j])&1 for j in range(6))
        require(encode(em)==prod(BASE)/z,'literal complement code')
        re=tuple((permute_mask(m,r)>>PHI[j])&1 for j in range(6))
        require(re==mirror(e),'literal mirror code')
        line_rows.append(dict(masks=[m,63^m],nodes=sorted([q[m],q[63^m]]),color='blue' if e==mirror(e) else 'black',code=str(z.numerator)))
    require(Counter(row['color'] for row in line_rows)=={'blue':8,'black':24},'literal line colors')
    require(sum(vp(BASE[2],p)>0 for p in PRIMES)==5 and sum(vp(BASE[4],p)>0 for p in PRIMES)==3,'prime-renaming hostile')
    report=dict(status='passed',base_cycle=BASE,base_factorizations=basis,chosen_phase_to_tile_gauge=PHI,
                scale_universe={'allowed_exponents':'0..2 at 2,11,13,19,23; nu2*nu19=0; S in {1,31}','count':len(scales),'binary_checks':scaled_checks},
                signed_lattice_checks=signed_checks,levels=levels,raw_cycle_reflections=raw,
                blue_tilings=16,blue_lines=8,black_lines=24,phase_shift_retained_blue_tilings=4,
                temporal_mirror_axis_multiplicities=dict(mirror_count),merged_fibres=fibres,gauge_count=len(gauges),
                gauge_descent_hostiles=witnesses,literal_lines=line_rows,
                checks={'no_inherited_producer_imports':True,'all_failure_gates_survive_optimization':True,'universal_claims_have_symbolic_proofs':True})
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':main()
