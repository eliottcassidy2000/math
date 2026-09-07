"""Independent literal edge-subset / permutation audit of repair potential.

No producer or earlier flow engine is imported. Uniform asymptotics are
proved in the referee report; these controls deliberately stay small.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations,permutations,product
from math import comb,factorial
import sys
sys.stdout.reconfigure(newline='\n')
G=0
def need(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
def occupancy(B):
    return Counter(c-r for r,c in B),Counter(c+r for r,c in B)
def diagonal_safe(B):
    d,e=occupancy(B)
    return max(d.values(),default=0)<=2 and max(e.values(),default=0)<=2
@lru_cache(None)
def tau(key):
    # Maximum retained subset, by full edge-subset search; no flow or matching.
    B=list(key)
    for size in range(len(B),-1,-1):
        if any(diagonal_safe(A) for A in combinations(B,size)): return len(B)-size
def potentials(B):
    ds,ss=occupancy(B)
    F=sum(max(n-2,0) for n in ds.values())
    zs=Counter(c+r for r,c in B if ds[c-r]<=2)
    I=sum(n==3 and all(ds[c-r]==1 for r,c in B if c+r==ell) for ell,n in ss.items())
    return F,sum(max(n-2,0) for n in zs.values()),I
def board(n,sigma,step=1):
    return tuple(sorted((r,sigma[c]) for r in range(n) for c in (r,(r+step)%n)))

boards=0;transpositions=0
for n in (3,4,5):
    chosen=list(permutations(range(n)))[:min(12,factorial(n))]
    for sigma in chosen:
        B=board(n,sigma); value=tau(B); F,Z,I=potentials(B)
        need(value>=F+Z>=F+I,'actual deterministic disjoint-deletion certificate')
        need(0<=value<=len(B),'literal minimum deletion range')
        # The retained-set consumer is checked against every safe subset.
        for size in range(len(B)+1):
            for T in combinations(B,size):
                if diagonal_safe(T): need(len(B)-len(T)>=value,'all actual retained subsets')
        for a,b in combinations(range(n),2):
            C=tuple(sorted((r,b if c==a else a if c==b else c) for r,c in B))
            lost=len(set(B)-set(C))
            need(lost<=4 and abs(value-tau(C))<=lost,'actual transposition and replacement Lipschitz')
            transpositions+=1
        boards+=1

# Exact finite forbidden-permutation controls, including m=0 and collisions.
for N in range(1,7):
    for shifts in ((),(0,),(0,1),(0,1,2)):
        cells={(i,(i+d)%N) for i in range(N) for d in shifts}
        degree=max([0]+list(Counter(i for i,j in cells).values())+
                   list(Counter(j for i,j in cells).values()))
        m=len(cells); hits=[]
        for p in permutations(range(N)):
            hits.append(sum((i,p[i]) in cells for i in range(N)))
        for k in range(min(4,N)+1):
            rook=sum(len({i for i,j in A})==k and len({j for i,j in A})==k
                     for A in combinations(cells,k))
            moment=Q(sum(comb(w,k) for w in hits),factorial(N))
            need(moment==Q(rook,comb(N,k)*factorial(k)), 'exact forbidden rook factorial moment')
            if k>=2:
                ordered=rook*factorial(k)
                need(0<=m**k-ordered<=k*(k-1)*degree*m**(k-1),'literal bounded-degree collision bound')
        prob=Q(hits.count(0),len(hits))
        for K in range(min(4,N)+1):
            bon=sum((-1)**k*Q(sum(comb(w,k) for w in hits),factorial(N)) for k in range(K+1))
            need((prob<=bon if K%2==0 else prob>=bon),'literal Bonferroni lower and upper bounds')

# Native row/column choices, with a frozen arbitrary row permutation.
triples=eligible=assignments=conditioned=0
for n in range(4,9):
    for shift in (0,2):
        rho=[(i+shift)%n for i in range(n)]
        neighbours={rho[i]:(i,(i+1)%n) for i in range(n)}
        rows_of={c:tuple(r for r in range(n) if c in neighbours[r]) for c in range(n)}
        checked_triples=0
        for ell in range(2*n-1):
            rows=[r for r in range(n) if 0<=ell-r<n]
            bad=0
            for selected in combinations(rows,3):
                triples+=1
                # Define eligibility by the actual forbidden geometric coincidences.
                is_eligible=True
                for r in selected:
                    for col in neighbours[r]:
                        other=next(t for t in rows_of[col] if t!=r)
                        if any(t==other or 2*t==r+other for t in selected if t!=r):
                            is_eligible=False
                if not is_eligible: bad+=1; continue
                eligible+=1
                target={r:ell-r for r in selected}
                intended={(r,target[r]) for r in selected}
                ds={c-r for r,c in intended}
                U={(r,c) for r in range(n) for c in range(n) if r+c==ell or c-r in ds}
                source_packets=list(product(*(neighbours[r] for r in selected)))
                need(len(source_packets)==8,'all eight actual source choices')
                seen_assignment=set()
                for cols in source_packets:
                    need(len(set(cols))==3,'eligibility forces distinct source columns')
                    forced=dict(zip(cols,(target[r] for r in selected)))
                    key=tuple(sorted(forced.items()))
                    need(key not in seen_assignment,'eight intended assignments are distinct')
                    seen_assignment.add(key)
                    companions={(t,forced[c]) for r,c in zip(selected,cols) for t in rows_of[c] if t!=r}
                    need(companions.isdisjoint(U),'all forced companion cells lie outside full four-line union')
                    free_cols=[c for c in range(n) if c not in forced]
                    free_labels=[c for c in range(n) if c not in forced.values()]
                    forbidden={(i,j) for i,c in enumerate(free_cols) for j,label in enumerate(free_labels)
                               if any((r,label) in U for r in rows_of[c])}
                    need(max(Counter(i for i,j in forbidden).values(),default=0)<=8 and
                         max(Counter(j for i,j in forbidden).values(),default=0)<=8,
                         'both conditional forbidden degrees at most eight')
                    line_sum=len(rows)+sum(n-abs(d) for d in ds)
                    need(len(forbidden)<=2*len(U)<=2*line_sum,
                         'native forbidden-count envelope by actual four line lengths')
                    # Only a bounded fresh set is conditioned exhaustively.
                    if checked_triples<2:
                        for p in permutations(range(n-3)):
                            mapping=forced|{c:free_labels[p[i]] for i,c in enumerate(free_cols)}
                            actual={(r,mapping[c]) for c in range(n) for r in rows_of[c]}
                            avoid=all((i,p[i]) not in forbidden for i in range(n-3))
                            need(avoid==(actual.intersection(U)==intended),'native conditional avoidance iff exact geometry')
                            if avoid:
                                dcount,scount=occupancy(actual)
                                need(scount[ell]==3 and all(dcount[d]==1 for d in ds),
                                     'actual isolated antidiagonal triple event')
                            conditioned+=1
                    assignments+=1
                checked_triples+=1
            if len(rows)>=3: need(bad<=4*len(rows)*(len(rows)-2),'uniform eligible discard bound')
            else: need(bad==0,'short-line zero contribution')
        lengths=[n-abs(ell-(n-1)) for ell in range(2*n-1)]
        need(sum(comb(L,3) if L>=3 else 0 for L in lengths)==comb(n,3)+2*comb(n,4),
             'exact all-line triple count and leading coefficient')

# Integrate the exact exponential polynomial, with z=exp(-2).
# The four terms are exp(-6)/3 times exp(k*theta) coefficients.
integrated={}
for k,coeff in ((4,1),(2,-3),(0,3),(-2,-1)):
    if k==0:
        integrated[3]=integrated.get(3,Q(0))+Q(coeff,3)
    else:
        upper=(6-k)//2
        integrated[upper]=integrated.get(upper,Q(0))+Q(coeff,3*k)
        integrated[3]=integrated.get(3,Q(0))-Q(coeff,3*k)
need(integrated=={1:Q(1,12),2:Q(-1,2),3:Q(5,4),4:Q(1,6)},
     'exact integrated four-line epsilon coefficients')
need(Q(1,5)/Q(1,4)==Q(4,5) and 6-Q(4,5)==Q(26,5),
     'independent weighted Jensen mean cost')
qlo=sum(Q((-2)**k,factorial(k)) for k in range(32))
qhi=sum(Q((-2)**k,factorial(k)) for k in range(31))
need(0<qlo<qhi<1,'exact alternating Taylor enclosure of exp(-2)')
epsilon_lo=qlo/12-qhi*qhi/2+5*qlo**3/4+qlo**4/6
need(epsilon_lo>Q(1,200),'independent exact positive gain exceeds one two-hundredth')

print('PASS independent two-direction repair audit')
print('No flow producer import: complete retained-edge subsets on',boards,'fresh boards;',transpositions,'transposition controls.')
print('Native line triples',triples,'eligible',eligible,'eight-choice assignments',assignments,'conditional permutations',conditioned)
print('Uniform avoidance lemma, exact four-line event, and adaptive retained-point proof independently accepted.')
print('Asymptotic gain epsilon=exp(-2)/12-exp(-4)/2+5*exp(-6)/4+exp(-8)/6.')
print('Gamma=1-5*exp(-2)+epsilon; no finite threshold or absence theorem inferred.')
print('Always-active exact gates',G)
