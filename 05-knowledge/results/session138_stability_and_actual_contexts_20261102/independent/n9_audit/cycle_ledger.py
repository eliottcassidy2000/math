"""Exact symbolic controls for the order-nine analytic reduction.
No producer imports; all failures remain active under python -O.
"""
from collections import defaultdict
from fractions import Fraction as Q
import json

def need(b, message):
    if not b:
        raise RuntimeError(message)

def parts(n, lo=2):
    yield ()
    for j in range(lo, n+1):
        for tail in parts(n-j, j):
            yield (j,)+tail

def clean(d):
    return {k:v for k,v in d.items() if v}

def add(*polys):
    d=defaultdict(Q)
    for p in polys:
        for k,v in p.items(): d[k]+=v
    return clean(d)

def scale(p,c): return clean({k:c*v for k,v in p.items()})
def literal(spec):
    return {tuple(int(i) for i in k):Q(v) for k,v in spec.items()}

def generic(n):
    universe=list(parts(n))
    need(len(universe)==len(set(universe)), 'partition uniqueness')
    H={t:Q(2**len(t)) for t in universe if all(j%2 for j in t)}
    u={t:Q((-1)**sum(j-1 for j in t)) for t in universe}
    v={t:Q((-1)**len(t)) for t in universe if sum(t)==n}
    return H,u,v,universe

H,u,v,universe=generic(9)
HH=literal({'':1,'3':2,'5':2,'7':2,'9':2,'33':4,'35':4,'333':8})
uu=literal({'':1,'2':-1,'3':1,'4':-1,'5':1,'6':-1,'7':1,'8':-1,'9':1,
 '22':1,'23':-1,'24':1,'25':-1,'26':1,'27':-1,'33':1,'34':-1,'35':1,'36':-1,
 '44':1,'45':-1,'222':-1,'223':1,'224':-1,'225':1,'233':-1,'234':1,'333':1,
 '2222':1,'2223':-1})
vv=literal({'9':-1,'27':1,'36':1,'45':1,'225':-1,'234':-1,'333':-1,'2223':1})
need(H==HH and u==uu and v==vv, 'complete cycle ledgers')
X=literal({'2':1,'22':-1,'222':1,'2222':-1,'3':1,'4':1,'5':1,'6':1,'7':1,'8':1,
 '23':1,'24':-1,'25':1,'26':-1,'33':3,'34':1,'35':3,'44':-1,'223':-1,
 '224':1,'233':1})
Gm=add(H,scale(u,-1),scale(v,-1)); Gp=add(H,scale(u,-1),v)
need(Gm==add(X,literal({'9':2,'333':8})), 'Gminus signs')
need(Gp==add(X,literal({'27':2,'36':2,'45':2,'225':-2,'234':-2,'333':6,'2223':2})), 'Gplus signs')

def select_cycle(j,local):
    d=defaultdict(Q)
    for t,c in local.items():
        global_type=tuple(sorted((j,)+t))
        d[global_type]+=global_type.count(j)*c
    return clean(d)

_,u5,v5,_=generic(5)
Cm=select_cycle(4,add(u5,v5,literal({'4':Q(1,2)})))
Cp=select_cycle(4,add(u5,scale(v5,-1),literal({'4':Q(1,2)})))
need(Cm==literal({'4':1,'24':-1,'34':1,'44':-1,'224':1}), 'C4 minus multiplicities')
need(Cp==add(Cm,literal({'45':2,'234':-2})), 'C4 plus multiplicities')
T=select_cycle(3,literal({'':1,'3':Q(3,2),'22':-1}))
need(T==literal({'3':1,'33':3,'223':-1}), 'triangle complement multiplicity')
C6=literal({'6':1,'26':-1}); C5p=literal({'5':1,'225':-2}); C5m=literal({'5':1})
P=literal({'2':1,'22':-1,'222':1,'2222':-1})
keep=add(P,literal({'23':1,'7':1,'8':1}))
remaining={}
for name,G,C4,C5 in [('minus',Gm,Cm,C5m),('plus',Gp,Cp,C5p)]:
    rem=add(G,scale(add(C4,T,C6,C5,keep),-1))
    need(all(c>=0 for c in rem.values()), 'unused products nonnegative: '+name)
    remaining[name]={''.join(map(str,k)):str(v) for k,v in rem.items()}
# Exact completion of the square, read as coefficients [constant,V,V^2].
quad=[Q(0),Q(5,8),Q(-1,3)]
square=[Q(75,256)-Q(15,16)**2/3,Q(2,3)*Q(15,16),Q(-1,3)]
need(quad==square,'six-vertex square completion')
need(1-Q(75,256)==Q(181,256),'triangle complement reserve')
need(Q(181,256)-Q(19,64)==Q(105,256),'final triangle reserve')
need(1-2*3*Q(1,16)==Q(5,8),'five-cycle reserve')
need(1-3*Q(1,4)==Q(1,4),'six-cycle reserve')
# Endpoint coefficient scaling, not an enumeration of graph inequalities.
scaled=[1024*(1-Q(175,256))/4,-Q(1024,16),Q(1024,64),-Q(1024,256),
        1024*Q(19,64)/16,Q(1024,64)]
need(scaled==[81,-64,16,-4,19,16], 'graph endpoint scale')
witness=[5,9,7,2,1,3]
need(sum(c*x for c,x in zip(scaled,witness))==0,'P3 plus 3K2 witness')
need(sum(c*x for c,x in zip(scaled,witness))-19==-19,'wedge necessity')
need(1-Q(3,4)**4==Q(175,256),'four fair pair sharpness')
print(json.dumps({'status':'PASS','cycle_types_including_empty':len(universe),
 'signed_ledgers':2,'cycle_complement_groups':3,'unused_nonnegative_products':remaining,
 'endpoint_residual_coefficients':[int(c) for c in scaled],
 'sharp_graph_witness_counts':witness,'variance_coefficient':'175/256',
 'triangle_reserve':'105/256'},sort_keys=True,indent=2))
