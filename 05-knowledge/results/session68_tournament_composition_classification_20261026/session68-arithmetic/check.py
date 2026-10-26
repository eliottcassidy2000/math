#!/usr/bin/env python3
"""Bounded independent arithmetic/composition checks. Standard library only.
No repository scripts are imported or executed. See REPORT.md for all-order proofs.
"""
from itertools import product, combinations, permutations
from collections import Counter
import json, hashlib

def require(x,msg):
    if not x: raise RuntimeError(msg)

def carousel(m):
    return [[int(0<(j-i)%m<=m//2) for j in range(m)] for i in range(m)]

def modules(a):
    n=len(a); out=[]
    for mask in range(1,1<<n):
        vs=[v for v in range(n) if mask>>v&1]
        if all(all(a[x][y]==a[x][vs[0]] for y in vs) for x in range(n) if not mask>>x&1): out.append(mask)
    return out

def aut_brute(a):
    n=len(a)
    return sum(all(a[i][j]==a[p[i]][p[j]] for i in range(n) for j in range(i)) for p in permutations(range(n)))

def expand(word):
    m=len(word); q=carousel(m); a=[[0]*4 for _ in range(4)]
    for j in range(1,4): a[0][j]=1
    a[1][2]=a[2][3]=a[3][1]=1
    return [[q[i//4][j//4] if i//4!=j//4 else (a[j%4][i%4] if word[i//4] else a[i%4][j%4]) for j in range(4*m)] for i in range(4*m)]

def triangles(a):
    return frozenset((i,j,k) for i,j,k in combinations(range(len(a)),3) if a[i][j]+a[j][k]+a[k][i] in (0,3))

def hamilton(a):
    n=len(a); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if dp[mask][v]:
                for w in range(n):
                    if not mask>>w&1 and a[v][w]: dp[mask|1<<w][w]+=dp[mask][v]
    return sum(dp[-1])

def aut_backtrack(a):
    # Independently uses actual expanded adjacency, without supplied block partition.
    n=len(a); deg=list(map(sum,a)); candidates=[[j for j in range(n) if deg[j]==deg[i]] for i in range(n)]
    order=sorted(range(n),key=lambda i:(len(candidates[i]),i)); assigned={}; used=set()
    def rec(k):
        if k==n:return 1
        v=order[k]; total=0
        for w in candidates[v]:
            if w in used or any(a[v][u]!=a[w][z] for u,z in assigned.items()):continue
            assigned[v]=w;used.add(w);total+=rec(k+1);used.remove(w);del assigned[v]
        return total
    return rec(0)

def period(w):
    n=len(w)
    return next(d for d in range(1,n+1) if n%d==0 and all(w[i]==w[(i+d)%n] for i in range(n)))

def canonical(w):return min(w[i:]+w[:i] for i in range(len(w)))

def mu(n):
    f=0;p=2
    while p*p<=n:
        if n%p==0:
            n//=p;f+=1
            if n%p==0:return 0
            while n%p==0:n//=p
        p+=1
    return (-1)**(f+(n>1))

def primitive(d):return sum(mu(e)*2**(d//e) for e in range(1,d+1) if d%e==0)

receipt={'pin':'fc3768ac5fcd87022fde9e6d7e3f06c0e5bc22e8','quotients':[], 'necklaces':[], 'anchor':[], 'crt':{}}
for m in (3,5,7,9,11):
    q=carousel(m); mods=modules(q)
    require(len(mods)==m+1,'carousel not prime')
    for i in range(m):
        counts={j:sum(q[i][z]*q[j][z] for z in range(m)) for j in range(m) if q[i][j]}
        require(max(counts,key=counts.get)==(i+1)%m,'successor not intrinsic')
    ac=aut_brute(q) if m<=7 else None
    require(ac is None or ac==m,'full quotient automorphism mismatch')
    receipt['quotients'].append({'m':m,'nonempty_modules':len(mods),'full_aut_brute':ac})
    words=list(product((0,1),repeat=m)); counts=Counter(period(w) for w in {canonical(w) for w in words})
    require(all(counts[d]==primitive(d)//d for d in counts),'Mobius class count failed')
    require(sum(d*c for d,c in counts.items())==2**m,'orbit mass failed')
    # Direct labelled C3-fibre verification for EVERY word through m=7.
    tri=triangles(expand(words[0]))
    require(len(tri)==(8*m**3-5*m)//3,'triangle count failed')
    if m<=7:
        require(all(triangles(expand(w))==tri for w in words),'labelled C3 set changed')
    receipt['necklaces'].append({'m':m,'vertices':4*m,'labelled_words':2**m,'isomorphism_classes':sum(counts.values()),'strata':[{'least_period':d,'classes':counts[d],'aut_order':3**m*m//d} for d in sorted(counts)],'labelled_C3_count':len(tri),'all_labelled_C3_sets_checked':m<=7})
for w in product((0,1),repeat=3):
    a=expand(w); ac=aut_backtrack(a); h=hamilton(a)
    require(ac==3**3*3//period(w),'expanded Aut mismatch')
    require(h==462915,'anchor Hamilton count mismatch')
    mods=modules(a); strong=[s for s in mods if all(not(s&t and s&~t and t&~s) for t in mods)]
    maximal=[s for s in strong if s!=(1<<12)-1 and not any(s!=t and s&t==s and t!=(1<<12)-1 for t in strong)]
    require(set(maximal)=={15<<(4*i) for i in range(3)},'blocks not canonical')
    receipt['anchor'].append({'word':list(w),'aut':ac,'H':h,'maximal_strong_modules':maximal})
# CRT test copied as literal witness, not via an inherited script.
w=(0,9,1,16,15,14,29,21,30)
lines=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
local={p:tuple(x%p for x in w) for p in (5,7,35)}
collision={p:{(i,j) for i,j in combinations(range(9),2) if local[p][i]==local[p][j]} for p in local}
for p in local:
    require(all(x in {r*r%p for r in range(p)} for x in local[p]),'not all square residues')
    require(len({sum(local[p][i] for i in line)%p for line in lines})==1,'not modular magic')
require(collision[35]==collision[5]&collision[7]==set(),'CRT collision identity failed')
require(len(set(w))==9,'CRT entries not distinct')
require(len({sum(w[i] for i in line) for line in lines})>1,'integer magic hostile vanished')
receipt['crt']={'witness':w,'local_square_residue_counts':{p:len({r*r%p for r in range(p)}) for p in (5,7)},'collision_counts':{p:len(collision[p]) for p in local},'ordinary_line_sums':[sum(w[i] for i in line) for line in lines]}
# Odd reflection followed by bit complement has no fixed necklace.
receipt['converse'] = []
for m in (3,5,7,9,11):
    reps={canonical(w) for w in product((0,1),repeat=m)}
    dual=lambda w: canonical(tuple(1-w[-i%m] for i in range(m)))
    require(all(dual(w)!=w and dual(dual(w))==w and period(dual(w))==period(w) for w in reps),'complement-reflection hostile')
    receipt['converse'].append({'m':m,'ordinary_classes':len(reps),'merged_classes':len(reps)//2,'self_converse_classes':0})
# Explicit normal-basis F8 map. f=t^3+t+1, alpha=1+t.
def mul8(a,b):
    r=0
    while b:
        if b&1:r^=a
        b>>=1;a<<=1
        if a&8:a^=11
    return r
sq=lambda a:mul8(a,a)
basis=(3,sq(3),sq(sq(3)))
field=lambda w:__import__('functools').reduce(int.__xor__,(b for z,b in zip(w,basis) if z),0)
require(basis==(3,5,7) and len({field(w) for w in product((0,1),repeat=3)})==8,'F8 not normal basis')
rows=[]
for w in product((0,1),repeat=3):
    x=field(w);require(field(w[-1:]+w[:-1])==sq(x),'Frobenius equivariance failed')
    y=sq(x);d=1
    while y!=x:y=sq(y);d+=1
    require(d==period(w),'field degree mismatch')
    trace=x^sq(x)^sq(sq(x));norm=mul8(mul8(x,sq(x)),sq(sq(x)))
    rows.append({'word':list(w),'field_polynomial_bits':x,'degree':d,'trace':trace,'norm':norm,'aut_order':81//d})
require(rows[4]['trace']==rows[7]['trace']==rows[4]['norm']==rows[7]['norm']==1 and rows[4]['aut_order']!=rows[7]['aut_order'],'trace/norm hostile failed')
bad_basis=(2,sq(2),sq(sq(2)))
require(bad_basis==(2,4,6) and bad_basis[0]^bad_basis[1]^bad_basis[2]==0,'bad-basis hostile absent')
reflection=lambda x: next(field(tuple(w[-i%3] for i in range(3))) for w in product((0,1),repeat=3) if field(w)==x)
nonmultiplicative=next((a,b) for a in range(8) for b in range(8) if reflection(mul8(a,b))!=mul8(reflection(a),reflection(b)))
receipt['finite_field']={'modulus_bits':11,'normal_basis_bits':basis,'rows':rows,'bad_basis_bits':bad_basis,'reflection_nonmultiplicative_witness':nonmultiplicative}
# Four-color tournament palette B_(e,f)=C3[A^e,A^f,1].
def palette(e,f):
    blocks=(expand((e,)),expand((f,)),[[0]])
    labels=[(i,j) for i,b in enumerate(blocks) for j in range(len(b))]
    return [[blocks[i][u][v] if i==j else int((j-i)%3==1) for j,v in labels] for i,u in labels]
pals={(e,f):palette(e,f) for e in (0,1) for f in (0,1)}
# A distinct singleton fixes all rotations of the C3 quotient, so colors cannot interchange.
base_tri=triangles(pals[(0,0)])
receipt['palette']=[]
for (e,f),a in pals.items():
    mods=modules(a);strong=[z for z in mods if all(not(z&t and z&~t and t&~z) for t in mods)]
    maximal=[z for z in strong if z!=511 and not any(z!=t and z&t==z and t!=511 for t in strong)]
    require(set(maximal)=={15,240,256},'palette canonical modules failed')
    ac=aut_backtrack(a);require(ac==9,'palette aut failed')
    require(triangles(a)==base_tri and len(base_tri)==18,'palette labelled C3 fibre failed')
    receipt['palette'].append({'color':[e,f],'aut':ac,'cyclic_triples':len(base_tri),'maximal_strong_modules':maximal})
# Relative-phase law independently tested on explicit cycles, not theorem formula.
receipt['phase_gluing']=[]
for m in (3,5,9,15):
    divisors=[d for d in range(1,m+1) if m%d==0]
    for d in divisors:
        for e in divisors:
            # A single 1 in a primitive period block, then repeat.
            a=tuple(int(i%d==0) for i in range(m));b=tuple(int(i%e==0) for i in range(m))
            oa={a[i:]+a[:i] for i in range(m)};ob={b[i:]+b[:i] for i in range(m)}
            pairs={canonical(tuple(zip(x,y))) for x in oa for y in ob}
            from math import gcd,lcm
            require(len(pairs)==gcd(d,e),'relative phase count mismatch')
            require(all(period(w)==lcm(d,e) for w in pairs),'joint period mismatch')
            receipt['phase_gluing'].append({'m':m,'d':d,'e':e,'joint_orbits':len(pairs),'joint_period':lcm(d,e)})
print(json.dumps(receipt,indent=2,sort_keys=True))
