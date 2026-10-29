#!/usr/bin/env python3
"""Independent finite exact audit; no producer imports, no assert-dependent guards."""
import json
from fractions import Fraction as F
from itertools import product, combinations, permutations
from functools import lru_cache
COUNTS={}
def need(ok,kind):
    COUNTS[kind]=COUNTS.get(kind,0)+1
    if not ok: raise RuntimeError(kind)
def natural(n,positive=False):
    if type(n) is not int or n < int(positive): raise ValueError('integer domain')
def native(b,n):
    natural(b,True); natural(n,True)
    if b<3 or b%2!=1 or n%2!=1: raise ValueError('odd native domain')
    z=b*n+1;a=0
    while z%2==0: z//=2;a+=1
    return z,a
@lru_cache(None)
def energy_rank(n,q):
    natural(n,True); natural(q)
    if n%2!=1: raise ValueError('odd source')
    if q==0:return 0
    p=(n+1)*F(3,2)**q;k=0
    while p>=3: p*=F(27,28);k+=1
    return 3*k+q

def word_data(b,word):
    natural(b,True)
    if b<3 or b%2!=1:raise ValueError('odd multiplier')
    c=d=0
    for a in word:
        natural(a,True); c=b*c+2**d;d+=a
    return b**len(word),c,d

def tags(V,delta):
    if not V or len(set(V))!=len(V):raise ValueError('finite nonempty distinct speeds')
    for v in V:natural(v,True)
    if not 0<delta<F(1,2):raise ValueError('delta')
    return {(v,k,e):((F(k)+e*delta)/v)%1 for v in V for k in range(v) for e in (-1,1)}
def safe(t,V,d):
    return all(d <= (t*v)%1 <= 1-d for v in V)
def accepted(V,d):return {tag:t for tag,t in tags(V,d).items() if safe(t,V,d)}

# One-step rank checks include every leaf, branch, and direct-terminal edge
# for the stated rectangle, independently of an eventual-run search.
for n in range(1,1024,2):
    out,a=native(3,n)
    for q in range(1,7):
        qq=q+(2 if a%3==0 else -1)
        need(energy_rank(out,qq)<energy_rank(n,q),'native_rank_edges')
        ratio=F(out+1,n+1)*F(3,2)**(qq-q)
        need(ratio <= (F(27,28) if a%3==0 else 1),'native_potential_edges')
need([energy_rank(n,q) for n,q in [(13,1),(5,3),(1,2),(1,1),(1,0)]]==[163,162,38,4,0],'displayed_rank')
need(energy_rank(1,1)>energy_rank(1,0),'reset_is_not_decrease')
for n0 in range(1,128,2):
    for q0 in range(1,5):
        n,q,j=n0,q0,0;bound=energy_rank(n,q)
        while q:
            n,a=native(3,n);q+=2 if a%3==0 else -1;j+=1
            need(j<=bound,'full_return_bound')
for q0 in range(1,21):
    n,q=83,q0
    for i in range(3):n,a=native(5,n);q+=2 if a%5==0 else -1;need(q>0,'five_cycle_active')
    need((n,q)==(83,q0),'five_cycle_closure')

# Exact carry cylinder: replay two representatives of every short word;
# neighboring odd residues must fail exact replay, not merely divisibility.
for b in (3,5,7):
    for length in range(1,5):
        for word in product(range(1,5),repeat=length):
            slope,c,d=word_data(b,word);mod=2**(d+1)
            residue=((2**d-c)*pow(slope,-1,mod))%mod
            for n0 in (residue,residue+mod):
                n=n0;cc=dd=jj=0;L=0
                for a in word:
                    out,actual=native(b,n);need(a==actual,'exact_cylinder_replay')
                    Lnew=b*L+2**dd*((2**a-b)*n0-1)
                    cc=b*cc+2**dd;dd+=a;jj+=1;n=out
                    need(2**dd*n==b**jj*n0+cc,'carry_identity')
                    need(Lnew==2**dd*(n0-n),'debt_recurrence');L=Lnew
                    need((L>0)==(n<n0),'original_source_payment')
            hostile=(residue+2)%mod;n=hostile;actual=[]
            for _ in word:n,a=native(b,n);actual.append(a)
            need(tuple(actual)!=word,'cylinder_neighbor_rejected')
need(word_data(5,(1,1,5))==(125,39,7),'five_affine_hostile')
need(word_data(3,(1,1,5))==(27,19,7),'three_rational_hostile')

# All orientations of K4, independently test positive circulation via
# directed simple cycles and positive functional via all topological orders.
edges=list(combinations(range(4),2));dag=0
for mask in range(64):
    arcs=[e if mask>>i&1 else e[::-1] for i,e in enumerate(edges)]
    cyclic=False
    for k in (3,4):
        for cyc in permutations(range(4),k):
            cycle_edges=list(zip(cyc,cyc[1:]+cyc[:1]))
            if all(e in arcs for e in cycle_edges):
                cyclic=True;boundary=[0]*4
                for u,v in cycle_edges:boundary[u]-=1;boundary[v]+=1
                need(boundary==[0]*4,'positive_cycle_kernel')
    orders=[p for p in permutations(range(4)) if all(p.index(u)<p.index(v) for u,v in arcs)]
    need(bool(orders)==(not cyclic),'graph_alternative')
    if orders:
        dag+=1;f={v:i for i,v in enumerate(orders[0])}
        for u,v in arcs:need(f[v]-f[u]>=1,'positive_functional')
need(dag==24,'K4_transitive_count')

# Endpoint existence is checked independently by partitioning the full circle
# into endpoints and open cells; every cell midpoint has constant feasibility.
birth_increase=False
for size in range(1,6):
    for V in combinations(range(1,7),size):
        for d in (F(1,3),F(1,4),F(1,7)):
            acc=accepted(V,d);bounds=sorted(set(tags(V,d).values()))
            midpoints=[((x+(bounds[(i+1)%len(bounds)]+(i==len(bounds)-1)))/2)%1 for i,x in enumerate(bounds)]
            feasible=any(safe(t,V,d) for t in bounds+midpoints)
            need(bool(acc)==feasible,'endpoint_existence')
            for w in set(range(1,7))-set(V):
                final=accepted(V+(w,),d)
                old={tag:t for tag,t in acc.items() if safe(t,(w,),d)}
                new={tag:t for tag,t in tags((w,),d).items() if safe(t,V,d)}
                need(not (old.keys() & new.keys()),'tag_disjointness')
                need(final==old|new,'birth_filter')
                birth_increase |= len(final)>len(acc)
need(birth_increase,'accepted_count_can_increase')
acc=accepted(tuple(range(1,14)),F(1,14))
need(len(acc)==12 and len(set(acc.values()))==6,'LRC14_isolated_fixture')
need(len(tags(tuple(range(1,14)),F(1,14)))==182,'LRC14_candidates')

# Finite aggregate budget and sharp strictness checks (not evidence of an
# infinite Collatz budget). Truncated costs equal summed survival at each K.
weights=[F(1,2),F(1,3),F(1,6)];times=[0,3,11]
M=sum(w*t for w,t in zip(weights,times))
for K in range(20):
    A=sum(w for w,t in zip(weights,times) if t>K)
    trunc=sum(w*min(t,K) for w,t in zip(weights,times))
    need(trunc==sum(sum(w for w,t in zip(weights,times) if t>s) for s in range(K)),'budget_survival_identity')
    for w,t in zip(weights,times):
        need(not A<w or t<=K,'atomic_strict_extraction')
        need(t<=M//w,'aggregate_time_bound')
need(sum(w for w,t in zip(weights,times) if t>10)==weights[2] and times[2]>10,'equality_not_sufficient')
# Independently test bounded-reset compiler on hostile nonmonotone budgets.
# The source coordinate is unchanged on local edges and strictly drops at resets.
for B in ([0,0,0,0,0,0], [0,9,0,8,1,7], [0,1,4,9,16,25]):
    prefix=[sum(B[m]+1 for m in range(n)) for n in range(len(B))]
    for n in range(len(B)):
        for r in range(B[n]+1):
            for rr in range(r):need(prefix[n]+rr<prefix[n]+r,'compiled_local_edges')
            for nn in range(n):
                for rr in range(B[nn]+1):need(prefix[nn]+rr<prefix[n]+r,'compiled_paid_resets')
need(native(3,3)==(5,1),'return_without_source_descent')
# malformed-input checks remain enabled with python -O
for fn,args in [(native,(3,0)),(native,(4,3)),(native,(3,2)),(native,(3,True)),(energy_rank,(3,-1)),(word_data,(3,(0,))), (tags,((),F(1,4))),(tags,((0,),F(1,4))),(tags,((1,1),F(1,4))),(tags,((1,),F(1,2)))]:
    try: fn(*args)
    except ValueError:COUNTS['malformed_rejected']=COUNTS.get('malformed_rejected',0)+1
    else:raise RuntimeError('malformed input accepted')
print(json.dumps({'status':'PASS','counts':COUNTS,'uses_assert':False,'producer_imports':False,'scope':'finite exact audit; infinite claims require written proof'},indent=2,sort_keys=True))
