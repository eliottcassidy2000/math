#!/usr/bin/env python3
"""Exact, standard-library-only checks; no assert gates or external inputs."""
from fractions import Fraction
from itertools import product
from math import prod, gcd
import json

P = (2,3,5,7,11,13,17,19,23,29)
F = ((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
LETTERS = 'ABCDEFGHIJKLMN'
BLOCK = (1,3,2,3,2,2,3,2,1,0,0,3,1,1)
PHI = (-4,0,-4,2,2,2,2,3,-1,3)
TOKEN = (1,0,1,0,0,0,0,0,1,0)
checks=0

def need(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)

def val(n,p):
    a=0
    while n%p == 0:
        a+=1
        n//=p
    return a

NUM=tuple(tuple(val(a,p) for p in P) for a,b in F)
DEN=tuple(tuple(val(b,p) for p in P) for a,b in F)
DELTA=tuple(tuple(a-b for a,b in zip(u,v)) for u,v in zip(NUM,DEN))

def dot(a,b): return sum(x*y for x,y in zip(a,b))
def number(x): return prod(p**e for p,e in zip(P,x))
def choose(x):
    for k,d in enumerate(DEN):
        if all(a>=b for a,b in zip(x,d)): return k
    raise RuntimeError('Missing denominator-one instruction')
def step(x):
    k=choose(x)
    return k,tuple(a+b for a,b in zip(x,DELTA[k]))
def integer_step(n):
    for k,(a,b) in enumerate(F):
        if n%b == 0: return k,n//b*a
    raise RuntimeError('Missing denominator-one instruction')
def rank(rows):
    a=[list(map(Fraction,row)) for row in rows]
    r=0
    for c in range(len(a[0])):
        j=next((j for j in range(r,len(a)) if a[j][c]),None)
        if j is None: continue
        a[r],a[j]=a[j],a[r]
        q=a[r][c]
        a[r]=[v/q for v in a[r]]
        for j in range(len(a)):
            if j!=r and a[j][c]:
                q=a[j][c]
                a[j]=[v-q*w for v,w in zip(a[j],a[r])]
        r+=1
    return r

def controller(a,b,c,d,s):
    return (a,b,c,d)+tuple(int(p==s) for p in P[4:])

for k in range(14):
    need(gcd(*F[k])==1, 'Fractions must be reduced')
    need(max(NUM[k]+DEN[k])<=1, 'Binary consumption/production premise')
    need(dot(PHI,DELTA[k])==(0,-2,1,0)[BLOCK[k]], 'Height grading')
    need(dot(TOKEN,DELTA[k])==int(k==13), 'Monotone token')
need(rank(DELTA)==10,'Full increment rank')
# Unknowns (ten character weights, q for block1, r for block2).
constraints=[list(v)+[-int(BLOCK[k]==1),-int(BLOCK[k]==2)] for k,v in enumerate(DELTA)]
need(rank(constraints)==11,'One-dimensional grading solution space')
need(all(dot(row,PHI+(-2,1))==0 for row in constraints),'Primitive grading solution')
need(gcd(*PHI)==1,'Primitive character')
need([k for k,v in enumerate(DELTA) if dot(PHI,v)==0 and not any(v[:4])]==[9,10], 'Horizontal pure-control pair')

# Four free nonnegative integer parameters g=D,j=J,t=C,m=M.
def circulation(g,j,t,m):
    return (g+2*t+m,g+t+m,t,g,g+t+2*m,g+t+2*m,g,t,0,j,t+j,t+m,m,0)
generators=[circulation(1,0,0,0),circulation(0,1,0,0),circulation(0,0,1,0),circulation(0,0,0,1)]
need(rank(generators)==4,'Four independent circulation generators')
for g,j,t,m in product(range(4),repeat=4):
    z=circulation(g,j,t,m)
    need(all(x>=0 for x in z),'Nonnegative circulation')
    need(all(sum(z[k]*DELTA[k][r] for k in range(14))==0 for r in range(10)), 'Exact circulation')
    counts=[sum(z[k] for k in range(14) if BLOCK[k]==b) for b in range(4)]
    need(counts[2]==2*counts[1], 'Closed-word height balance')
    need(counts[3]-counts[1]==2*g,'Closed-word second count balance')

# All first-two-step behaviors: cap every exponent at 2. No sampling.
graph={b:set() for b in range(4)}
witness={}
for x in product(range(3),repeat=10):
    k,y=step(x)
    ell,z=step(y)
    edge=(BLOCK[k],BLOCK[ell])
    graph[edge[0]].add(edge[1])
    n=number(x)
    if edge not in witness or n<witness[edge]['integer']:
        witness[edge]={'integer':n,'instructions':LETTERS[k]+LETTERS[ell]}
    ik,ny=integer_step(n)
    il,nz=integer_step(ny)
    need((ik,il,ny,nz)==(k,ell,number(y),number(z)),'Independent integer/vector check')
need(graph=={0:{0,1,2},1:{0,1,2,3},2:{0,1,2,3},3:{0,1,2,3}},'Complete global successor graph')

onehot={b:set() for b in range(4)}
onehot_witness={}
for a,b,c,d in product(range(3),repeat=4):
    for s in (1,11,13,17,19,23,29):
        x=controller(a,b,c,d,s)
        k,y=step(x)
        ell,z=step(y)
        edge=(BLOCK[k],BLOCK[ell])
        onehot[edge[0]].add(edge[1])
        n=number(x)
        if edge not in onehot_witness or n<onehot_witness[edge]['integer']:
            onehot_witness[edge]={'integer':n,'instructions':LETTERS[k]+LETTERS[ell]}
        need(sum(y[4:])<=1 and sum(z[4:])<=1,'One-hot invariance')
need(onehot=={0:{0,1,2},1:{0,1,2,3},2:{0,2,3},3:{0,1,2,3}},'Complete one-hot successor graph')
need((choose(controller(0,1,0,1,11)),choose(controller(0,1,0,1,13)))==(4,0),'231/273 phase witness')
need(integer_step(247)==(7,1001) and integer_step(1001)==(0,187),'Multimarker forbidden-edge stress test')

# Complete primitive-six-cycle census. For each ABDEFG permutation, each
# coordinate is consumed at most once, so cap2 preserves every prefix guard.
six_words={}
g_state_cases={'no17_forcesH':0,'positive7_forcesC':0,'multiple17_forcesC':0,'six_cycle':0}
for x in product(range(3),repeat=10):
    y=x
    word=''
    for _ in range(6):
        k,y=step(y)
        word+=LETTERS[k]
    if y==x and set(word)==set('ABDEFG'):
        six_words[word]=six_words.get(word,0)+1
    predicted=(x[1]==x[2]==x[3]==x[9]==0 and x[6]==1 and x[4]>=1 and x[8]>=1 and x[0]*x[7]==0)
    need((word=='GBDEFA')==predicted,'Exact GBDEFA guard family')
    if word[0]=='G' and x[5]+x[6]>0 and x[4]+x[9]>0 and x[7]+x[8]>0:
        if x[6]==0:
            q=x
            for _ in range(2*x[0]+x[8]+2):
                k,q=step(q)
                if LETTERS[k] not in 'DG': break
            need(LETTERS[k]=='H','No17 branch must leave periodic subalphabet at H')
            g_state_cases['no17_forcesH']+=1
        elif x[3]>0:
            need(word.startswith('GBAC'),'Positive7 branch must fire C')
            g_state_cases['positive7_forcesC']+=1
        elif x[6]>=2:
            need(word.startswith('GBC'),'Multiple17 branch must fire C')
            g_state_cases['multiple17_forcesC']+=1
        else:
            need(word=='GBDEFA' and y==x,'Only periodic G-state branch')
            g_state_cases['six_cycle']+=1

rotations={'GBDEFA'[k:]+'GBDEFA'[:k] for k in range(6)}
need(set(six_words)==rotations,'Only one cyclic six-instruction ordering')
base_cycle=[4301,17765,16302,9867,8671,23023]
need(gcd(*base_cycle)==1,'Primitive scaling base')
for n,k,nn in zip(base_cycle,'GBDEFA',base_cycle[1:]+base_cycle[:1]):
    need(integer_step(n)==(LETTERS.index(k),nn),'Actual minimal three-control cycle')
for a,r,v,w,s in product(range(3),repeat=5):
    scale=2**a*11**r*13**v*19**w*23**s
    legal=all(integer_step(n*scale)==(LETTERS.index(k),nn*scale) for n,k,nn in zip(base_cycle,'GBDEFA',base_cycle[1:]+base_cycle[:1]))
    need(legal==(a*w==0),'Scaled six-cycle family and hostile products')

# Weak first-exit observation: remove J,K, retain silent divergence.
def weak_exit(x):
    seen=set()
    while x not in seen:
        seen.add(x)
        k,y=step(x)
        if BLOCK[k]!=0: return BLOCK[k],y
        x=y
    return None
for a,b,c,d in product(range(4),repeat=4):
    left=weak_exit(controller(a,b,c,d,11))
    right=weak_exit(controller(a,b,c,d,13))
    need((left==right)==(b*d==0),'Exact guarded reciprocal contraction')
    if b*d:
        need((left[0],right[0])==(2,1),'Interior phase separation')
    elif b==d==0:
        need(left is None,'Silent divergence')

for a in range(51):
    x=(a,)+(0,)*9
    for _ in range(a):
        k,x=step(x)
        need(k==11,'Power-of-two countdown')
    need(choose(x)==13,'First exit after exact countdown')

# Full-order trace sanity, not evidence substituted for all-time proofs.
x=(1,)+(0,)*9
prime_events=[]
for t in range(37012):
    if t and all(e==0 for e in x[1:]): prime_events.append([t,x[0]])
    _,x=step(x)
need(prime_events==[[19,2],[69,3],[281,5],[710,7],[2375,11],[3893,13],[8102,17],[11361,19],[19268,23],[36981,29]],'Inherited first-ten output control')

out={
 'status':'PASS',
 'explicit_checks':checks,
 'fraction_order':[f'{a}/{b}' for a,b in F],
 'blocks':{str(b):''.join(LETTERS[k] for k in range(14) if BLOCK[k]==b) for b in range(4)},
 'height_weights':dict(zip(map(str,P),PHI)),
 'height_changes':[0,-2,1,0],
 'grading_constraint_rank':11,
 'increment_rank':10,
 'nonnegative_circulation_generators':{name:dict((LETTERS[k],v) for k,v in enumerate(z) if v) for name,z in zip(('g','j','t','m'),generators)},
 'global_successor_graph':{str(b):sorted(v) for b,v in graph.items()},
 'global_witnesses':{f'{a}->{b}':witness[(a,b)] for a,b in sorted(witness)},
 'onehot_successor_graph':{str(b):sorted(v) for b,v in onehot.items()},
 'onehot_witnesses':{f'{a}->{b}':onehot_witness[(a,b)] for a,b in sorted(onehot_witness)},
 'complete_threshold_universe':3**10,
 'complete_onehot_threshold_universe':7*3**4,
 'primitive_six_cycle_words_capped_counts':six_words,
 'global_periodic_G_state_elimination_cases':g_state_cases,
 'proved_global_period_set':[2,6],
 'primitive_six_cycle_base':base_cycle,
 'primitive_six_cycle_scale':'2^a 11^r 13^v 19^w 23^s times spectator factor; a*w=0',
 'guarded_merge_witness':[231,273],
 'multimarker_edge_witness':[247,1001,187],
 'countdown_family_checked':[0,50],
 'seed2_prime_events':prime_events,
 'scope':'Finite checks verify complete two-step threshold universes; general grading, kernel, guarded merge and infinite-memory proofs are in README.md.'
}
print(json.dumps(out,sort_keys=True,indent=2))
