#!/usr/bin/env python3
"""Independent finite controls; no producer imports, no optimization-disabled asserts."""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, gcd
import json

checks = 0
def require(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(label)

def rank(rows, p=0):
    if not rows:
        return 0
    a = [[Fraction(x) if not p else x % p for x in row] for row in rows]
    r = 0
    for j in range(len(a[0])):
        z = next((z for z in range(r, len(a)) if a[z][j]), None)
        if z is None:
            continue
        a[r], a[z] = a[z], a[r]
        inv = 1 / a[r][j] if not p else pow(a[r][j], -1, p)
        a[r] = [x * inv if not p else x * inv % p for x in a[r]]
        for z in range(r + 1, len(a)):
            v = a[z][j]
            if v:
                a[z] = [x - v*y if not p else (x-v*y) % p for x,y in zip(a[z],a[r])]
        r += 1
        if r == len(a):
            break
    return r

def compose(a,b):
    return tuple(a[x] for x in b)

def closure(n, generators):
    identity=tuple(range(n)); seen={identity}; todo=[identity]
    while todo:
        a=todo.pop()
        for b in generators:
            c=compose(a,b)
            if c not in seen:
                seen.add(c); todo.append(c)
    return sorted(seen)

def sign(seq):
    return (-1)**sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))

def exterior_action(g, s):
    values=[g[x] for x in s]
    return tuple(sorted(values)), sign(values)

def check_monomial(basis, group, generators, action, label):
    index={s:i for i,s in enumerate(basis)}
    unseen=set(basis); good=bad=0
    while unseen:
        s=min(unseen)
        orbit={action(g,s)[0] for g in group}
        negative=any(action(g,s)==(s,-1) for g in group)
        bad+=negative; good+=not negative
        unseen-=orbit
    rows=[]
    for g in generators:
        for s in basis:
            target, eps=action(g,s)
            row=[0]*len(basis); row[index[target]]+=eps; row[index[s]]-=1
            if any(row): rows.append(row)
    for p in (0,2,3,5):
        actual=len(basis)-rank(rows,p)
        require(actual==good+(bad if p==2 else 0), (label,p,actual,good,bad))
    if len(group)%2:
        require(bad==0, (label,'odd-group sign obstruction'))
    return {'label':label,'dimension':len(basis),'group_order':len(group),'oriented_orbits':good,'nonoriented_orbits':bad}

records=[]
for n in range(1,8):
    cycle=tuple((i+1)%n for i in range(n))
    gens=[cycle]
    g=closure(n,gens)
    for r in range(n+1):
        records.append(check_monomial(list(combinations(range(n),r)),g,gens,exterior_action,f'C{n}/exterior{r}'))
for n in range(2,6):
    cycle=tuple((i+1)%n for i in range(n)); trans=list(range(n));trans[0],trans[1]=1,0
    gens=[cycle,tuple(trans)]; g=closure(n,gens)
    for r in range(n+1):
        records.append(check_monomial(list(combinations(range(n),r)),g,gens,exterior_action,f'S{n}/exterior{r}'))
gens=[tuple(3*((i//3+1)%3)+i%3 for i in range(9)),tuple(3*(i//3)+(i%3+1)%3 for i in range(9))]
g=closure(9,gens)
for r in range(10):
    records.append(check_monomial(list(combinations(range(9),r)),g,gens,exterior_action,f'C3xC3/exterior{r}'))

# A separate signed action: endpoint reversal, NOT edge-subset permutation sign.
for n in range(2,5):
    edges=list(combinations(range(n),2)); edge_index={e:i for i,e in enumerate(edges)}
    basis=list(range(1<<len(edges)))
    def graph_action(g,mask):
        target=0; reversals=0
        for j,(a,b) in enumerate(edges):
            if mask>>j&1:
                aa,bb=g[a],g[b]
                reversals+=aa>bb
                target|=1<<edge_index[tuple(sorted((aa,bb)))]
        return target,(-1)**reversals
    cycle=tuple((i+1)%n for i in range(n));trans=list(range(n));trans[0],trans[1]=1,0
    gens=[cycle,tuple(trans)]; group=list(permutations(range(n)))
    records.append(check_monomial(basis,group,gens,graph_action,f'S{n}/graph-Fourier'))

# Sign-convention hostile: swapping endpoints of a single fixed edge has
# graph orientation sign -1 but induced permutation of its one edge is +1.
require(exterior_action((0,), (0,))[1]==1,'single-edge induced permutation sign')
require(sign([1,0])==-1,'endpoint reversal sign')

# Walsh integral-lattice certificate. F^T F=N I; modulo2 all columns agree.
walsh=[]
for d in range(1,7):
    n=1<<d
    f=[[(-1)**((a&b).bit_count()) for b in range(n)] for a in range(n)]
    for a in range(n):
        for b in range(n):
            require(sum(f[z][a]*f[z][b] for z in range(n))==(n if a==b else 0),('Walsh',d,a,b))
    require(rank(f,2)==1,('Walsh rank2',d))
    require(rank(f,3)==n,('Walsh rank3',d))
    walsh.append({'d':d,'size':n,'index_power_of_2':d*(1<<(d-1)),'rank_mod2':1})

# General base-change failure outside the signed permutation class:
# order3 action on the rank2 augmentation ideal.
aug=[[-1,-1],[1,-2]] # g-I for g=[[0,-1],[1,-1]].
require(rank(aug)==2 and rank(aug,3)==1 and rank(aug,2)==2,'augmentation ideal defect at3')

# Integer jet certificate of the 25-gap example (not tournament realization).
p=[(comb(25,i)+(24*comb(5,i//5) if i%5==0 else 0))//25 for i in range(26)]
for i in range(26):
    require(25*p[i]==comb(25,i)+(24*comb(5,i//5) if i%5==0 else 0),('integer coefficient',i))
jets=[sum((-1)**(i-j)*comb(i,j)*p[i] for i in range(j,26)) for j in range(26)]
depth=lambda prime:next(j for j,c in enumerate(jets) if c and (not prime or c%prime))
require([depth(a) for a in (0,5,2)]==[5,21,25],'three cancellation depths')
require(all((p[i]-comb(25,i))%2==0 for i in range(26)),'odd-action mod2 binomial')
def valuation(n,p):
    if n==0:return None
    v=0
    while n%p==0:n//=p;v+=1
    return v
jet_record={'coefficients':p,'jets_at_minus_one':jets,'v5_jets':[valuation(x,5) for x in jets],'depth_Q':depth(0),'depth_mod5':depth(5),'depth_mod2':depth(2)}
gcd_chain=[]; content=0
for c in jets:
    content=gcd(content,c); gcd_chain.append(content)
require(gcd_chain==[0]*5+[3000]*8+[600]*4+[120]*4+[24]*4+[1],'exact jet gcd filtration')
jet_record['gcd_filtration_Delta_1_to_26']=gcd_chain
jet_record['depth_mod_5_power']=[next(j for j,c in enumerate(jets) if c%(5**e)) for e in range(1,5)]
require(jet_record['depth_mod_5_power']==[21,17,13,5],'jet prime-power filtration')

# C3: direct exterior multiplication, not a Hilbert-series inference.
def wedge(a,b):
    result={}
    for s,cs in a.items():
        for t,ct in b.items():
            if set(s)&set(t):continue
            key=tuple(sorted(s+t))
            result[key]=result.get(key,0)+cs*ct*sign(s+t)
    return {s:c for s,c in result.items() if c}
x={(0,):1,(1,):1,(2,):1}
y={(0,1):1,(1,2):1,(0,2):-1}
require(wedge(x,y)=={(0,1,2):3},'C3 xy=3z')
require(wedge(y,x)=={(0,1,2):3},'C3 yx=3z')
require(wedge(x,x)=={} and wedge(y,y)=={},'C3 other products')
require(rank([[3]],3)==0 and rank([[3]],2)==1,'C3 product characteristic boundary')
small_jets=[sum((-1)**(i-j)*comb(i,j) for i in range(j,4)) for j in range(4)]
require(small_jets==[0,2,-2,1],'C3 scalar jets')
require([12-6+1,-6+2,1]==[7,-4,1],'q=1 translated torification coefficients')

# Smith exponents are invisible over fields: [2] and [4] lose one rank mod2,
# but the annihilator counts over Z/4 are2 and4.
require(rank([[2]],2)==rank([[4]],2)==0,'same field rank')
require(sum(2*x%4==0 for x in range(4))==2,'annihilator2 overZ4')
require(sum(4*x%4==0 for x in range(4))==4,'annihilator4 overZ4')

print(json.dumps({'status':'PASS','checks':checks,'monomial_cases':len(records),'records':records,'walsh':walsh,'clone_jet_control':jet_record,'C3_product_coefficient':3,'C3_jets':small_jets,'controls':['C3 augmentation ideal has genuine mod3 invariant defect','C3 exterior invariants commute but multiplication changes mod3 and scalar depth changes mod2','single-edge sign conventions differ','Walsh lattice not integral at2','Smith exponent invisible to field ranks','negative torification coefficient is not torsion']},indent=2))
