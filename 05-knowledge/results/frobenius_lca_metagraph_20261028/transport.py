"""Relative Frobenius phase: exact joint tables, query costs, mixed periods."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
import json
from math import gcd, lcm, prod
from pathlib import Path
from fields import Field, need, rank


def query_costs(matrices):
    r=len(matrices)
    queries=sorted({sum(row[j]<<k for k,row in enumerate(matrices)) for j in range(r*r)})
    @lru_cache(None)
    def cost(states):
        if states.bit_count()<=1: return 0
        values=[]
        for query in queries:
            yes=states&query; no=states^yes
            if yes and no: values.append(1+max(cost(yes),cost(no)))
        need(bool(values),'phase queries separate states')
        return min(values)
    domain=(1<<r)-1
    nonadaptive=next(k for k in range(r+1) if any(len({tuple((q>>i)&1 for q in qs) for i in range(r)})==r for qs in combinations(queries,k)))
    return cost(domain),nonadaptive


def phase_queries(data):
    r=data['degree']; f=Field(r,data['polynomial_bits']); q=1<<(r//2)
    normal=data['normal_basis']
    def B(a,b): return f.trace(f.mul(a,f.power(b,q)))
    matrices=[]
    for k in range(r):
        matrices.append(tuple(B(a,f.power(b,1<<k)) for a in normal for b in normal))
    need(len(set(matrices))==r,'faithful relative phase table')
    need(all(B(f.mul(a,a),f.mul(b,b))==B(a,b) for a in normal for b in normal),'common phase gauge')
    adaptive,nonadaptive=query_costs(matrices)
    # A linear readout of a whole joint row is a stronger query than one entry.
    # Invertible circulant rows can encode arbitrary binary labels of phases.
    decoders=[]
    for bit in range((r-1).bit_length()):
        weights=next(w for w in range(1<<r) if all(
            sum((w>>j&1)*matrices[k][j] for j in range(r))%2==((k>>bit)&1)
            for k in range(r)))
        decoders.append(weights)
    for k in range(r):
        reconstructed=sum((sum((w>>j&1)*matrices[k][j] for j in range(r))%2)<<bit
                          for bit,w in enumerate(decoders))
        need(reconstructed==k,'joint parity phase decoder')
    spectrum={}; generators=0; representatives=[]
    for a in range(1,f.size):
        orbit=f.orbit(a)
        if len(orbit)!=r or rank(orbit)!=r: continue
        generators+=1
        if a==min(orbit): representatives.append(a)
        row=tuple(B(a,b) for b in orbit)
        if row not in spectrum:
            shifted=[tuple(row[(j+k-i)%r] for i in range(r) for j in range(r)) for k in range(r)]
            ca,cn=query_costs(shifted)
            spectrum[row]=dict(trace_row=row,adaptive=ca,nonadaptive=cn,normal_generators=0)
        spectrum[row]['normal_generators']+=1
    pair_costs=Counter(); adaptive_pair_costs=Counter(); aggregate_pairs=0
    for a,b in product(representatives,repeat=2):
        row=tuple(B(a,f.power(b,1<<j)) for j in range(r))
        shifted=[tuple(row[(j+k-i)%r] for i in range(r) for j in range(r)) for k in range(r)]
        ca,cn=query_costs(shifted); adaptive_pair_costs[ca]+=1; pair_costs[cn]+=1
        # Different normal orbit pairs still have invertible joint row tables.
        need(rank([sum(x<<j for j,x in enumerate(v[:r])) for v in shifted])==r,'cross-normal invertibility')
        for bit in range((r-1).bit_length()):
            need(any(all(sum((w>>j&1)*shifted[k][j] for j in range(r))%2==((k>>bit)&1)
                         for k in range(r)) for w in range(1<<r)),'cross-orbit parity decoder')
        aggregate_pairs+=1
    if r==6: need(dict(pair_costs)=={3:8,5:8},'incoming F64 orbit-pair census')
    if r==8: need(dict(pair_costs)=={4:192,6:32,7:32},'incoming F256 orbit-pair census')
    full_period_hostile=None
    if r in (6,8):
        a,b=2,(34 if r==6 else 18); c=b^1
        profile=lambda y:sum(B(f.power(a,1<<j),y)<<j for j in range(r))
        need(len(f.orbit(a))==r and len(f.orbit(b))==r and f.power(b,1<<(r//2))==c,'full-period hostile phase')
        need(profile(b)==profile(c) and rank(f.orbit(a))==r-1,'normality guard hostile')
        full_period_hostile=dict(anchor=a,targets=[b,c],common_profile=profile(b),anchor_span_rank=r-1)
    return dict(degree=r,trace_row=matrices[0][:r],distinct_relative_tables=len(set(matrices)),
                binary_storage_bits=(r-1).bit_length(),minimum_adaptive_incidence_queries=adaptive,
                minimum_nonadaptive_incidence_queries=nonadaptive,
                optimal_joint_parity_queries=len(decoders),joint_parity_decoder_masks=decoders,
                all_normal_basis_generators=generators,normal_basis_query_spectrum=list(spectrum.values()),
                full_period_not_normal_hostile=full_period_hostile,
                different_orbit_scope=dict(normal_orbits=len(representatives),ordered_pairs=aggregate_pairs,
                    nonadaptive_cost_histogram=dict(sorted(pair_costs.items())),
                    adaptive_cost_histogram=dict(sorted(adaptive_pair_costs.items())),
                    aggregate_query_cost=(r-1).bit_length()),
                all_joint_tables=matrices)


def mixed_periods(m):
    divisors=[d for d in range(1,m+1) if m%d==0]
    edges=list(combinations(range(3),2)); cases=0; graph_tests=0; worst_min=0
    examples=[]
    for ds in product(divisors,repeat=3):
        states=list(product(*(range(d) for d in ds)))
        L=lcm(*ds)
        def canon(x): return min(tuple((x[i]+s)%ds[i] for i in range(3)) for s in range(L))
        classes={x:canon(x) for x in states}
        need(len(set(classes.values()))==prod(ds)//L,'CRT joint class count')
        good=[]
        for mask in range(8):
            selected=[e for j,e in enumerate(edges) if mask>>j&1]
            fibres={}
            for x in states:
                sig=tuple((x[j]-x[i])%gcd(ds[i],ds[j]) for i,j in selected)
                fibres.setdefault(sig,set()).add(classes[x])
            separates=all(len(v)==1 for v in fibres.values())
            layers=[]
            for pa in range(2,m+1):
                primes=[p for p in range(2,pa+1) if pa%p==0 and all(p%d for d in range(2,p))]
                if len(primes)!=1: continue
                layer={i for i,d in enumerate(ds) if d%pa==0}
                if not layer: continue
                seen={min(layer)}
                changed=True
                while changed:
                    old=len(seen)
                    for i,j in selected:
                        if i in layer and j in layer and (i in seen or j in seen): seen.update((i,j))
                    changed=len(seen)!=old
                layers.append(seen==layer)
            need(separates==all(layers),'prime-layer theorem versus complete orbit enumeration')
            if separates: good.append(mask)
            graph_tests+=1
        minimum=min(x.bit_count() for x in good); worst_min=max(worst_min,minimum)
        if ds in [(2,3,6),(2,6,6),(2,4,8),(4,8,8)]:
            examples.append(dict(periods=ds,joint_classes=prod(ds)//L,minimum_edges=minimum,
                                 minimal_graphs=[x for x in good if x.bit_count()==minimum]))
        cases+=1
    # The short-period centre loses a phase that a maximum-period centre sees.
    ds=(2,m,m); x=(0,0,0); y=(0,0,2)
    star=lambda z:tuple((z[j]-z[0])%gcd(ds[0],ds[j]) for j in (1,2))
    need(star(x)==star(y) and (x[2]-x[1])%m!=(y[2]-y[1])%m,'short anchor hostile')
    return dict(field_degree=m,period_tuples=cases,graph_tests=graph_tests,largest_minimum_edges=worst_min,
                examples=examples,short_anchor_hostile=dict(periods=ds,states=[x,y],same_star_record=star(x)))


def main():
    p=argparse.ArgumentParser(); p.add_argument('--output'); args=p.parse_args()
    data=json.loads(Path(__file__).with_name('fields-normal.json').read_text())
    result=dict(relative_phase=[phase_queries(data[str(1<<r)]) for r in (4,6,8)],
                mixed_periods=[mixed_periods(r) for r in (6,8)])
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: Path(args.output).write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')


if __name__=='__main__': main()
