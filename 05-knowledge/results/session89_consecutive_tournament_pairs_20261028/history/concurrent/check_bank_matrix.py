#!/usr/bin/env python3
from itertools import permutations,product
from collections import Counter
from math import isqrt,prod
import json

checks=0

def require(ok,message):
    global checks
    checks+=1
    if not ok:raise RuntimeError(message)

prime_grid=((19,59,2),(7,83,41),(89,139,13))
edge_grid=((1,4,5),(34,2,32),(16,24,8))
basis=(58,20,1,26,6,11)
probes=(48,29,25,56)
L=(1,2,12,7)
edges=((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))

def gf_mul(a,b):
    ans=0
    while b:
        if b&1:ans^=a
        b>>=1;a<<=1
        if a&64:a^=67
    return ans

def gf_pow(x,k):
    ans=1
    while k:
        if k&1:ans=gf_mul(ans,x)
        x=gf_mul(x,x);k>>=1
    return ans

def trace(x):
    ans=0
    for _ in range(6):ans^=x;x=gf_mul(x,x)
    require(ans in (0,1),('trace',ans))
    return ans

def linear(x,cols):
    ans=0
    for k,c in enumerate(cols):
        if x>>k&1:ans^=c
    return ans

def source_probes(mask):
    deg=[0]*4
    for k,(i,j) in enumerate(edges):
        if mask>>k&1:deg[i]^=1;deg[j]^=1
    return sum(deg[i]<<i for i in range(3))|((mask.bit_count()%2)<<3)

require(len({linear(x,basis) for x in range(64)})==64,'field basis invertible')
require(len({linear(x,L) for x in range(16)})==16,'four map invertible')
for mask in range(64):
    x=linear(mask,basis)
    observed=sum(trace(gf_mul(x,gf_pow(p,8)))<<i for i,p in enumerate(probes))
    require(observed==source_probes(mask),('trace probe identity',mask))

for i,j in product(range(3),repeat=2):
    p=prime_grid[i][j]
    require(p>=2 and all(p%d for d in range(2,isqrt(p)+1)),('actual prime',p))
    q=p+1;sfmask=0
    for k,prime in enumerate((2,3,5,7)):
        e=0
        while q%prime==0:q//=prime;e^=1
        sfmask|=e<<k
    require(q==1,('prime plus one support',p))
    require(sfmask==linear(source_probes(edge_grid[i][j]),L),('bank dictionary',i,j,p))
    require(gf_pow(linear(edge_grid[i][j],basis),9)==1,('norm one',i,j))

cuts={sum((a[i]^b[j])<<(3*i+j) for i,j in product(range(3),repeat=2))
      for a in product(range(2),repeat=3) for b in product(range(2),repeat=3)}
require(len(cuts)==32,'cut dimension five')

def rect(word):
    bit=lambda i,j:(word>>(3*i+j))&1
    return tuple(bit(0,0)^bit(i,0)^bit(0,j)^bit(i,j) for i,j in product((1,2),repeat=2))

def natural(word,p):
    inv=[p.index(i) for i in range(3)];out=0
    for i in range(3):out|=((word>>(3*inv[i]+inv[i]))&1)<<(3*i+i)
    for i in range(3):
        for j in range(i+1,3):
            u,v=inv[i],inv[j];lo,hi=sorted((u,v));flip=int(u>v)
            out|=(((word>>(3*lo+hi))&1)^flip)<<(3*i+j)
            out|=(((word>>(3*hi+lo))&1)^flip)<<(3*j+i)
    return out

def ordinary(word,p):
    return sum(((word>>(3*i+j))&1)<<(3*p[i]+p[j]) for i,j in product(range(3),repeat=2))

weights=Counter();natural_failures=0
for word in range(512):
    sigma=prod(prime_grid[i][j]+1 for i,j in product(range(3),repeat=2) if word>>(3*i+j)&1)
    square=isqrt(sigma)**2==sigma
    require(square==(word in cuts)==(not any(rect(word))),('square/cut/rectangles',word))
    if word in cuts:
        weights[word.bit_count()]+=1
        a=[(word>>(3*i))&1 for i in range(3)]
        b=[((word>>j)&1)^(word&1) for j in range(3)]
        require(word==sum((a[i]^b[j])<<(3*i+j) for i,j in product(range(3),repeat=2)),('reconstruct',word))
    for p in permutations(range(3)):
        require((ordinary(word,p) in cuts)==(word in cuts),('ordinary closure',word,p))
        natural_failures+=((natural(word,p) in cuts)!=(word in cuts))

action_rows=[]
for p in permutations(range(3)):
    surviving=sum(natural(word,p) in cuts for word in cuts)
    witness=next((word for word in sorted(cuts) if natural(word,p) not in cuts),None)
    action_rows.append({'permutation':p,'accepted_to_accepted':surviving,
                        'preserves_cut_code':surviving==32,
                        'zero_image_is_cut':natural(0,p) in cuts,
                        'first_failing_codeword':witness})
require([tuple(r['permutation']) for r in action_rows if r['preserves_cut_code']]==[(0,1,2)],'only identity preserves code')
require(sorted(r['accepted_to_accepted'] for r in action_rows)==[0,16,16,16,16,32],'all natural intersection counts')

hostile=natural(0,(1,0,2))
require(hostile==(1<<1)|(1<<3),'natural hostile slots')
require(rect(hostile)==(0,1,1,0),'natural hostile rectangles')
require(60*8==480 and isqrt(480)**2!=480,'natural hostile arithmetic')
require(any(rect(1)),'rank one E00 fails rectangles')
require(weights=={0:1,3:6,4:9,5:9,6:6,9:1},'weight enumerator')

print(json.dumps({'status':'PASS','always_active_checks':checks,'source_commit':'b38ec1a1eaf7f4e31a7a38308e2764b5c7650d6d','prime_grid':prime_grid,'edge_grid':edge_grid,'all_subset_universe':512,'square_sigma_subsets':len(cuts),'weight_enumerator':dict(sorted(weights.items())),'natural_action_predicate_mismatches_among_512_times_6':natural_failures,'natural_action_rows':action_rows,'hostile':{'source':0,'permutation':[1,0,2],'target':hostile,'prime_subset':[59,7],'sigma':480,'rectangle_syndrome':rect(hostile)},'scope':'Independent reconstruction of all supplied field/probe/prime data; no imported producer code; no replay of all eight wedge charts'},indent=2))
