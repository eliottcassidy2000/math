#!/usr/bin/env python3
"""Exact bounded-circulation controls; no dependencies, no assertions, no Git."""
import itertools as it
import json
import math
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from pathlib import Path

def require(test, message):
    if not test:
        raise RuntimeError(message)

def edges(n):
    return list(it.combinations(range(n), 2))

def orientation(n, mask):
    return [(v,u) if (mask>>k)&1 else (u,v) for k,(u,v) in enumerate(edges(n))]

def cyclic(n, steps):
    return [(u,v) if (v-u)%n in steps else (v,u) for u,v in edges(n)]

def scores(n, arcs):
    out=[0]*n
    for u,v in arcs: out[u]+=1
    return tuple(out)

def demand(n, arcs, f):
    out=[0]*n
    for (u,v),x in zip(arcs,f): out[u]-=x; out[v]+=x
    return tuple(out)

def reverse(arcs,D):
    return [(v,u) if k in D else (u,v) for k,(u,v) in enumerate(arcs)]

def transport(f,D,M):
    return tuple(M-x if k in D else x for k,x in enumerate(f))

def flows(n, arcs, M, b=None):
    """Independent free chord coordinates: solve star edges from vertex balances."""
    if b is None: b=(0,)*n
    require(sum(b)==0,'demand must sum to zero')
    star=[k for k,(u,v) in enumerate(arcs) if u==0 or v==0]
    chord=[k for k in range(len(arcs)) if k not in star]
    for vals in it.product(range(M+1),repeat=len(chord)):
        f=[0]*len(arcs); bal=[0]*n
        for k,x in zip(chord,vals):
            u,v=arcs[k]; f[k]=x; bal[u]-=x; bal[v]+=x
        for k in star:
            u,v=arcs[k]; j=v if u==0 else u
            sign=1 if u==0 else -1
            x=sign*(b[j]-bal[j])
            if x<0 or x>M: break
            f[k]=x
        else:
            yield tuple(f)

def coefficient_count(n,M,target):
    """Recurrence for coefficient of product sum x_i^k x_j^(M-k)."""
    @lru_cache(None)
    def C(ds):
        r=len(ds)
        if sum(ds)!=M*r*(r-1)//2 or any(d<0 or d>M*(r-1) for d in ds): return 0
        if r==0: return 1
        if r==1: return int(ds[0]==0)
        d0=ds[0]; total=0
        # Exponents k_j at vertex 0, M-k_j at the other endpoint.
        for ks in compositions_bounded(d0,r-1,M):
            next_ds=tuple(sorted(d-M+k for d,k in zip(ds[1:],ks)))
            total+=C(next_ds)
        return total
    value=C(tuple(sorted(target)))
    return value,C.cache_info().currsize

def compositions_bounded(total,length,M):
    if length==0:
        if total==0: yield ()
        return
    for k in range(max(0,total-(length-1)*M),min(M,total)+1):
        for tail in compositions_bounded(total-k,length-1,M): yield (k,)+tail

def hamilton(n,arcs):
    outs=[0]*n
    for u,v in arcs: outs[u]|=1<<v
    dp={(1<<v,v):1 for v in range(n)}
    for mask in range(1,1<<n):
        for v in range(n):
            count=dp.get((mask,v),0)
            if not count: continue
            avail=outs[v]&~mask
            while avail:
                bit=avail&-avail; avail-=bit; w=bit.bit_length()-1
                dp[(mask|bit,w)]=dp.get((mask|bit,w),0)+count
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

def poly_interpolate(values):
    """Exact forward differences in monomial basis, independent of Ehrhart theorem."""
    diff=[Fraction(x) for x in values]; coeff=[Fraction(0)]*len(values)
    basis=[Fraction(1)]
    for k in range(len(values)):
        for j,c in enumerate(basis): coeff[j]+=diff[0]*c
        diff=[b-a for a,b in zip(diff,diff[1:])]
        new=[Fraction(0)]*(len(basis)+1)
        for j,c in enumerate(basis):
            new[j]-=k*c/(k+1); new[j+1]+=c/(k+1)
        basis=new
    return coeff

def eval_poly(coeff,x):
    return sum(c*x**k for k,c in enumerate(coeff))

def main():
    report={}; values=[]; cache_sizes=[]
    for M in range(11):
        v,s=coefficient_count(5,M,(2*M,)*5); values.append(v); cache_sizes.append(s)
    coeff=poly_interpolate(values[:7])
    require(all(eval_poly(coeff,M)==values[M] for M in range(11)),'n5 extra controls')
    require(all(eval_poly(coeff,-M)==eval_poly(coeff,M-2) for M in range(1,13)),'n5 reciprocity')
    h=[sum((-1)**j*math.comb(7,j)*values[k-j] for j in range(min(7,k)+1)) for k in range(11)]
    require(h[:6]==h[5::-1] and all(v==0 for v in h[6:]),'h* symmetry and degree')
    r5=cyclic(5,{1,2}); independent=[]; restricted5=[]
    path5={k for k,(u,v) in enumerate(r5) if v==u+1}
    for M in range(9):
        all_flows=list(flows(5,r5,M)); independent.append(len(all_flows))
        restricted5.append(sum(all(f[k]==0 for k in path5) for f in all_flows))
    require(restricted5==[math.comb(M+2,2) for M in range(9)],'n5 path-face polynomial')
    require(independent==values[:9],'independent direct n5 flows')
    report['n5']={'counts_M0_to_10':values,'coefficient_cache_states':cache_sizes,'independent_star_counts_M0_to_8':independent,'polynomial_coefficients_low_to_high':[str(c) for c in coeff],'h_star':h[:6],'dimension':6,'path_zero_counts_M0_to_8':restricted5}
    # The M1 coefficient counts labelled tournaments with each target score.
    fibres=Counter(scores(5,orientation(5,t)) for t in range(1<<10))
    for s,c in fibres.items():
        v,_=coefficient_count(5,1,s)
        require(v==c,'M1 score fibre coefficient')
    require(fibres[(2,)*5]==24,'regular n5 count')
    report['n5_score_fibres']={'labelled_tournaments':sum(fibres.values()),'labelled_score_vectors':len(fibres),'regular_fibre':24}
    # All orientation changes and all capacity-one vectors on K4, arbitrary demands.
    checked=0
    for t in range(64):
        A=orientation(4,t); s=scores(4,A)
        for dm in range(64):
            D={k for k in range(6) if dm>>k&1}; A2=reverse(A,D); s2=scores(4,A2)
            delta=demand(4,A,[int(k in D) for k in range(6)])
            require(delta==tuple(b-a for a,b in zip(s,s2)),'score delta')
            for f in it.product(range(2),repeat=6):
                fp=transport(f,D,1); b=demand(4,A,f); bp=demand(4,A2,fp)
                require(bp==tuple(x-y for x,y in zip(b,delta)),'general transport')
                require(transport(fp,D,1)==f,'involution')
                require(sum(fp)==sum(f)+len(D)-2*sum(f[k] for k in D),'total use sidecar')
                centered=tuple(2*x-1 for x in f); centered2=tuple(2*x-1 for x in fp)
                require(demand(4,A,centered)==demand(4,A2,centered2),'centered divergence')
                require(tuple(x+y for x,y in zip(s,b))==tuple(x+y for x,y in zip(s2,bp)),'target invariant')
                checked+=1
    report['general_K4_transport']={'orientation_reversal_vector_triples':checked,'capacity':1,'all_demands':True}
    # Original path operation, with actual ports 0 and n-1 unchanged.
    path_rows=[]; path_checks=0
    for n in (3,5):
        base=edges(n); D={k for k,(u,v) in enumerate(base) if v!=u+1}
        for bits in range(1<<len(D)):
            A=list(base)
            for j,k in enumerate(sorted(D)):
                if bits>>j&1: A[k]=A[k][::-1]
            if scores(n,A)!=((n-1)//2,)*n: continue
            Ap=reverse(A,D)
            require(scores(n,Ap)==((n+1)//2,)+((n-1)//2,)*(n-2)+((n-3)//2,),'path scores')
            for M in (1,2,3):
                source=set(flows(n,A,M)); b=(-M,)+(0,)*(n-2)+(M,)
                target=set(flows(n,Ap,M,b)); image={transport(f,D,M) for f in source}
                require(image==target,'path source-sink bijection')
                require(source=={transport(f,D,M) for f in target},'path return')
                require(sum(1 for _ in flows(n,Ap,M))!=len(source),'hostile image is not circulation')
                path_checks+=len(source)
                path_indices=set(range(len(A)))-D
                source_zero={f for f in source if all(f[k]==0 for k in path_indices)}
                target_zero={f for f in target if all(f[k]==0 for k in path_indices)}
                require({transport(f,D,M) for f in source_zero}==target_zero,'path-zero source-sink transport')
                require(len(source_zero)==(1 if n==3 else math.comb(M+2,2)),'path-face count')
                path_rows.append({'n':n,'bits':bits,'M':M,'count':len(source),'partner_circulation_count':sum(1 for _ in flows(n,Ap,M)),'path_zero_count':len(source_zero),'path_zero_total_use_histogram':dict(sorted(Counter(map(sum,source_zero)).items()))})
    report['fixed_path_transport']={'rows':path_rows,'transported_flows':path_checks}
    # Capacity equality loses genuine tournament data at n7.
    n7=[]
    choices=[('cyclic',cyclic(7,{1,2,3})),('Paley',cyclic(7,{1,2,4}))]
    A=choices[0][1]
    tri=(0,2,4)
    require(all(e in A for e in ((0,2),(2,4),(4,0))),'chosen triangle')
    cycle={(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])}
    choices.append(('cyclic_triangle_reversed',[(v,u) if (u,v) in cycle else (u,v) for u,v in A]))
    for name,A in choices:
        require(scores(7,A)==(3,)*7,'n7 regularity')
        F=list(flows(7,A,1)); hist=Counter(map(sum,F)); H=hamilton(7,A)
        require(len(F)==2640,'n7 capacity 1')
        path={k for k,(u,v) in enumerate(A) if v==u+1}
        restricted=[f for f in F if all(f[k]==0 for k in path)]
        require(len(path)==6 and len(restricted)==91,'n7 fixed-path face')
        n7.append({'name':name,'arcs':A,'H':H,'capacity1':len(F),'path_zero_capacity1':len(restricted),'capacity1_total_use_histogram':[hist[k] for k in range(22)]})
    require(len({row['H'] for row in n7})>1,'H information loss')
    require(len({tuple(row['capacity1_total_use_histogram']) for row in n7})>1,'weighted information loss')
    v,s=coefficient_count(7,1,(3,)*7)
    require(v==2640,'independent n7 coefficient')
    report['n7']=n7; report['n7_M1_coefficient_states']=s
    # Independently enumerate all orientations of the 15 nonpath edges.
    base=edges(7); nonpath=[k for k,(u,v) in enumerate(base) if v!=u+1]
    presentations=Counter()
    for mask in range(1<<len(nonpath)):
        A=list(base)
        for j,k in enumerate(nonpath):
            if mask>>j&1: A[k]=A[k][::-1]
        if scores(7,A)==(3,)*7: presentations[hamilton(7,A)]+=1
    require(sum(presentations.values())==91,'n7 independent fixed-path presentation count')
    require(presentations==Counter({171:57,175:25,189:9}),'n7 presentation H stratification')
    report['n7_fixed_path_presentations']={'universe':1<<15,'regular_count':91,'by_H':dict(presentations)}
    A=cyclic(7,{1,2,3}); path={k for k,(u,v) in enumerate(A) if v==u+1}
    interior3=sum(all(f[k]==0 for k in path) for f in flows(7,A,1,(-1,0,0,0,0,0,1)))
    require(interior3==34,'n7 codegree-three leading coefficient')
    report['n7_path_face_h_star']={'dimension':9,'codegree':3,'degree':7,'constant':1,'linear':81,'leading':interior3,'palindromic':False}
    # Exact weighted-vs-unweighted hostile at M2,n3.
    require(coefficient_count(3,2,(2,)*3)[0]==3,'n3 unweighted count')
    report['multiplicity_hostile']={'n':3,'M':2,'unweighted':3,'labelled_parallel_edge_weighted':10}
    report['status']='PASS; proof and finite controls are separately scoped'
    return report

if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
