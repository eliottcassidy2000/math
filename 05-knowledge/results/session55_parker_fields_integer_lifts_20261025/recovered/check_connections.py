"""Exact independently implemented small controls; no external dependencies."""
import json
from math import isqrt
from itertools import product

def require(x, msg):
    if not x: raise RuntimeError(msg)
def lines(v):
    return [sum(v[3*r:3*r+3]) for r in range(3)]+[sum(v[c::3]) for c in range(3)]+[v[0]+v[4]+v[8],v[2]+v[4]+v[6]]
def normal(c,a,b):
    return [c+a,c-a-b,c+b,c-a+b,c,c+a-b,c-b,c+a+b,c-a]
def square(x): return x>=0 and isqrt(x)**2==x

def paths(n):
    adj={i:[j for j in range(1,n+1) if i!=j and square(i+j)] for i in range(1,n+1)}
    ans=[]
    def dfs(p,unused):
        if not unused:
            if p[0]<p[-1]: ans.append(p)
            return
        for v in adj[p[-1]]:
            if v in unused: dfs(p+[v],unused-{v})
    for s in range(1,n+1): dfs([s],set(range(1,n+1))-{s})
    return ans
out={}
# Universal linear normal form identities, with repeated-value hostile controls.
for c,a,b in product(range(-3,4), repeat=3):
    require(lines(normal(c,a,b))==[3*c]*8,'normal form')
    if a>b>0: require((len(set(normal(c,a,b)))==9)==(a!=2*b),'collision test')
out['normal_form_integer_controls']=7**3
# All binary-valued solutions, not just parametrized ones.
binary=[v for v in product(range(2),repeat=9) if len({s%2 for s in lines(v)})==1]
require(all(v[0]==v[8] and v[2]==v[6] and v[1]==v[3]==v[5]==v[7] for v in binary),'characteristic two duplication')
out['all_F2_magic_arrays']=len(binary)
# Four common-center APs do not by themselves give additive compatibility.
e=1105; c=e*e
D=[c-x*x for x in range(1,e) if square(2*c-x*x)]
completions=[(a,b) for a in D for b in D if a>b and a!=2*b and a-b in D and a+b in D]
require(len(D)>=4 and not completions,'common-center hostile')
out['common_center_hostile']={'root':e,'differences':sorted(D),'full_magic_pairs':completions}
# Independent plain DFS of transition counts, no endpoint forcing/pruning.
counts={}
for n in range(15,26):
    ps=paths(n); counts[n]=len(ps)
    if n==23:
        witness=next(p for p in ps if 22 in (p[0],p[-1]) and all(set(e)!={14,22} for e in zip(p,p[1:])))
        out['degree_two_endpoint_hostile_Q23']=witness
    if n==25:
        usage=sum(any(set(e)=={1,15} for e in zip(p,p[1:])) for p in ps)
        require(usage==0,'Q25 forbidden edge')
        out['Q25_edge_1_15_usage']=usage
require(list(counts.values())==[1,1,1,0,0,0,0,0,3,0,10],'15..25 census')
out['unoriented_path_counts']=counts
# Positive AP-square-sum inherited map; full graph positivity is a real filter.
for A,B,C in [(40,9,41),(63,16,65),(4,3,5)]:
    v=[2*(C*C-4*A*B),2*C*C,2*(C*C+4*A*B)]
    require(A*A+B*B==C*C and all(square(v[i]+v[j]) for i in range(3) for j in range(i)),'AP map')
out['AP_positive_control']=[482,3362,6242]
out['AP_positivity_hostile']=[-46,50,146]
print(json.dumps(out,indent=2,sort_keys=True))
