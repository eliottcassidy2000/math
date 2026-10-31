#!/usr/bin/env python3
"""Independent small merged-tournament census and deletion-coordinate audit.
Standard library; no inherited code or data imports. All failure checks survive -O.
Owner gauge: fixed path n->...->1, free bits in (x,y), y rising, x falling.
Canonical codes: upper triangular adjacency word in row-major order (MSB first).
"""
from itertools import permutations, combinations
from collections import Counter, defaultdict
from functools import lru_cache, reduce
from math import gcd, prod
import json

def require(c, message):
    if not c: raise RuntimeError(message)

@lru_cache(None)
def pairs(n): return tuple(combinations(range(n),2))
@lru_cache(None)
def perms(n): return tuple(permutations(range(n)))

def canonical(A):
    n=len(A); best=1<<len(pairs(n))
    for p in perms(n):
        c=0
        for i,j in pairs(n): c=(c<<1)|A[p[i]][p[j]]
        if c<best: best=c
    return best

def converse(A): return tuple(zip(*A))
def merged(A): return min(canonical(A), canonical(converse(A)))

def adjacency(n,mask):
    A=[[0]*n for _ in range(n)]
    for i in range(n-1): A[i][i+1]=1
    for k,(x,y) in enumerate(tiles(n)):
        i,j=n-x,n-y
        if (mask>>k)&1: A[j][i]=1
        else: A[i][j]=1
    return tuple(map(tuple,A))

def tiles(n): return tuple((x,y) for y in range(1,n-1) for x in range(n,y+1,-1))
def reflect(n,m):
    ts=tiles(n); ix={t:k for k,t in enumerate(ts)}
    return sum(((m>>k)&1)<<ix[(n-y+1,n-x+1)] for k,(x,y) in enumerate(ts))

def primitive(v):
    g=reduce(gcd,v)
    return tuple(x//g for x in v)

def census(n,previous=None):
    ts=tiles(n); full=(1<<len(ts))-1
    As=[adjacency(n,t) for t in range(full+1)]
    cs=[canonical(A) for A in As]
    ci={c:i for i,c in enumerate(sorted(set(cs)))}
    opposite={c:canonical(converse(A)) for c,A in zip(cs,As)}
    merges=sorted(set(min(c,opposite[c]) for c in cs))
    node_id={c:i for i,c in enumerate(merges)}
    ns=[node_id[min(c,opposite[c])] for c in cs]
    members=defaultdict(list)
    for m,u in enumerate(ns): members[u].append(m)
    nodes=[]
    for u,masks in sorted(members.items()):
        A=As[masks[0]]; p=perms(n)
        H=sum(all(A[q[i]][q[i+1]] for i in range(n-1)) for q in p)
        aut=sum(all(A[q[i]][q[j]]==A[i][j] for i,j in pairs(n)) for q in p)
        classes=sorted(set(ci[cs[m]] for m in masks))
        blue=[m for m in masks if reflect(n,m)==m]
        require(len(masks)*aut==H*len(classes),'Hamiltonian fibre law')
        row=dict(node=u,classes=classes,scores=sorted(map(sum,A)),H=H,aut=aut,
                 masks=masks,blue_masks=blue,type='pure_blue' if len(blue)==len(masks) else 'mixed' if blue else 'pure_black')
        if previous:
            counts=Counter(); counts2=Counter()
            for m in masks:
                B=As[m]
                L=tuple(tuple(B[i][j] for j in range(n-1)) for i in range(n-1))
                R=tuple(tuple(B[i][j] for j in range(1,n)) for i in range(1,n))
                counts[previous['node_codes'][merged(L)]]+=1
                counts2[previous['node_codes'][merged(R)]]+=1
            v=tuple(counts[k] for k in range(len(previous['nodes'])))
            require(counts==counts2,'Two endpoint-deletion profiles differ')
            row.update(face_vector=v,primitive_face_vector=primitive(v))
        nodes.append(row)
    lines=[]
    for m in range((full+1)//2):
        k=full^m
        lines.append(dict(masks=[m,k],nodes=sorted([ns[m],ns[k]]),color='blue' if reflect(n,m)==m else 'black'))
    edges=Counter((tuple(r['nodes']),r['color']) for r in lines)
    sc=sum(c==opposite[c] for c in ci)
    require(len(nodes)==(len(ci)+sc)//2,'Converse Burnside count')
    require(sum(len(x['masks']) for x in nodes)==1<<len(ts),'Cube partition')
    require(all(v%2==0 for (uv,c),v in edges.items() if c=='black'),'Black multiplicity parity')
    return dict(n=n,tiles=ts,ordinary_classes=len(ci),self_converse=sc,merged_nodes=len(nodes),
                nodes=nodes,lines=lines,edges=[dict(nodes=uv,color=c,multiplicity=v) for (uv,c),v in sorted(edges.items())],
                node_codes=node_id)

def main():
    data={}; prev=None
    for n in (3,4,5,6):
        now=census(n,prev)
        require((now['ordinary_classes'],now['self_converse'],now['merged_nodes'])=={3:(2,2,2),4:(4,2,3),5:(12,8,10),6:(56,12,34)}[n],'Inherited census mismatch')
        if prev:
            vecs=[tuple(r['primitive_face_vector']) for r in now['nodes']]
            require(len(set(vecs))==len(vecs),'Primitive face vectors not injective')
        data[str(n)]=now;prev=now
    # Independent labelled enumeration at n3..5 confirms no isomorphism class omitted by fixed-path atlas.
    for n in (3,4,5):
        cc=set()
        for b in range(1<<len(pairs(n))):
            A=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(pairs(n)):
                A[i][j]=(b>>k)&1; A[j][i]=1-A[i][j]
            cc.add(canonical(A))
        require(len(cc)==data[str(n)]['ordinary_classes'],'Independent full-labelled count mismatch')
    # Encode the 3-coordinate n5 fingerprints by arbitrary coordinate primes 2,3,5.
    for row in data['5']['nodes']:
        row['formal_prime_code_2_3_5']=prod(p**e for p,e in zip((2,3,5),row['primitive_face_vector']))
        a,b,c=row['primitive_face_vector']
        row['selected_code_7_17_55']=7**a*17**b*55**c
        row['selected_code_factorization']={'7':a,'17':b,'5':c,'11':c}
        z=row['selected_code_7_17_55']
        recovered=[]
        for p in (7,17,5):
            e=0
            while z%p==0:
                e+=1; z//=p
            recovered.append(e)
        require(tuple(recovered)==tuple(row['primitive_face_vector']), 'Prime valuation decoding failed')
        require(z==11**c,'Residual 11 valuation failed')
    require(len(set(row['selected_code_7_17_55'] for row in data['5']['nodes']))==10,'Selected code collision')
    # Minimal literal obstruction to strong lumpability of merged-node growth.
    lower_pair=data['4']['nodes'][1]['masks']
    sandwich=data['5']['nodes'][3]['masks']
    require(lower_pair==[2,4] and sandwich==[16],'Non-lumpability witness gauge changed')
    lift_counts=[sum(m>>3==a for m in sandwich) for a in lower_pair]
    require(lift_counts==[1,0],'Non-lumpability witness failed')
    data['non_lumpability_witness']={'lower_merged_node':1,'lower_masks':lower_pair,
        'upper_merged_node':3,'upper_masks':sandwich,'individual_lift_counts':lift_counts}
    axis_nodes=[r for r in data['6']['nodes'] if sum(x!=0 for x in r['primitive_face_vector'])==1]
    require([r['node'] for r in axis_nodes]==[0,5,17,32,33], 'Order6 axis census failed')
    require([r['node'] for r in axis_nodes if r['type']=='pure_blue']==[0,17], 'Order6 pure-blue axis census failed')
    S=adjacency(4,1)
    T=[[0]*6 for _ in range(6)]
    for j in range(1,6):T[0][j]=1
    for i in range(1,5):
        T[i][5]=1
        for j in range(1,5):T[i][j]=S[i-1][j-1]
    require(data['6']['node_codes'][merged(T)]==5,'Sandwich witness is not mixed axis node5')
    data['mixed_axis_counterexample']={'order':6,'node':5,'construction':'1 -> strong_order4 -> 1',
       'raw_profile':data['6']['nodes'][5]['face_vector'],
       'primitive_profile':data['6']['nodes'][5]['primitive_face_vector'],
       'fiber_size':len(data['6']['nodes'][5]['masks']),
       'blue_fiber_size':len(data['6']['nodes'][5]['blue_masks'])}
    # Do not publish a large literal n6 trace: its complete atlas is recomputed above;
    # this compact summary retains every tested primitive signature and the census.
    data['6'].pop('lines')
    data['6'].pop('edges')
    for row in data['6']['nodes']:
        row['fiber_size']=len(row.pop('masks'))
        row['blue_fiber_size']=len(row.pop('blue_masks'))
    data['checks']={'full_labelled_n3_n5':True,'primitive_profiles_injective_n4_n6':True,
      'coordinate_primes_are_external_labels':True,'failure_gates_survive_optimization':True}
    print(json.dumps(data,indent=2,sort_keys=True))
if __name__=='__main__':main()
