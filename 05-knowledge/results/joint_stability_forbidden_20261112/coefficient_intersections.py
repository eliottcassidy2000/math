"""Complete all-length overlap graph of the540 binary coefficient rays.

For different odd offsets B,D, any collision has
min(k,l)<=v2(B-D), because A*2^k-C*2^l=D-B.
Equal offsets give either disjoint rays or the same eventual dyadic ray.
"""
from itertools import combinations
from collections import defaultdict,Counter
import json
from census_nine import read,save,need

def power_exponent(x):return x.bit_length()-1 if x>0 and not x&(x-1) else None

def ray_count(A,B,X):
    q=(X-B)//A
    return max(0,q.bit_length()-1) if q>0 else 0

def run():
    laws=[tuple(r) for r in read('completed_bulk_heads.json.gz')['closed_laws']];classes=defaultdict(list)
    for i,(A,B) in enumerate(laws):classes[A//(A&-A),B].append(i)
    need(len(classes)==540,'no two declared laws are eventual-shift equivalent')
    coll=defaultdict(set);edges=[]
    for (i,(A,B)),(j,(C,D)) in combinations(enumerate(laws),2):
        hits=set()
        if B==D:continue
        d=abs(B-D);v=(d&-d).bit_length()-1
        for k in range(1,v+1):
            numerator=A*2**k+B-D
            if numerator%C==0:
                ell=power_exponent(numerator//C)
                if ell is not None and ell>=1:hits.add((k,ell,A*2**k+B))
            numerator=C*2**k+D-B
            if numerator%A==0:
                ell=power_exponent(numerator//A)
                if ell is not None and ell>=1:hits.add((ell,k,C*2**k+D))
        if hits:
            edges.append({'first':i,'second':j,'collisions':sorted(hits)})
            for k,ell,c in hits:coll[c].update([i,j])
    def count(X):return sum(ray_count(A,B,X) for A,B in laws)-sum(len(ids)-1 for c,ids in coll.items() if c<=X)
    bound=1<<20;literal=set()
    for A,B in laws:
        k=1
        while A*2**k+B<=bound:literal.add(A*2**k+B);k+=1
    need(count(bound)==len(literal),'independent literal image union')
    levels=[]
    for b in range(1,61):
        X=2**b;value=count(X)
        if b<=20:need(value==sum(c<=X for c in literal),'every finite prefix count')
        levels.append({'height_exponent':b,'distinct_coefficients_at_most_2^b':value})
    # Prove eventual dyadic-height linearity by monotone inequalities for
    # each individual ray, rather than extrapolating the finite table.
    thresholds=[];offset=0
    for A,B in laws:
        L=(A-1).bit_length();t=L+int(not A&(A-1) and B>0);b=t+1
        while not (B<=(2**t-A)*2**(b-t) and (2*A-2**t)*2**(b-t)>-B):b+=1
        need(2**t-A>=0 and 2*A-2**t>=0,'monotone height threshold inequalities')
        thresholds.append(b);offset+=t
    offset+=sum(len(ids)-1 for ids in coll.values())
    threshold=max(max(thresholds),max(coll).bit_length())
    need((threshold,offset)==(21,6078),'exact eventual linear count constants')
    need(all(count(2**b)==540*b-offset for b in range(threshold,61)),'proved-threshold finite controls')
    degrees=Counter()
    adjacency=[set() for _ in laws];pairvalues={}
    for e in edges:
        i,j=e['first'],e['second'];degrees[i]+=1;degrees[j]+=1;adjacency[i].add(j);adjacency[j].add(i);pairvalues[i,j]={c for k,l,c in e['collisions']}
    seen=set();components=[]
    for i in range(len(laws)):
        if i in seen:continue
        seen.add(i);stack=[i];component=[]
        while stack:
            v=stack.pop();component.append(v)
            for w in adjacency[v]-seen:seen.add(w);stack.append(w)
        components.append(sorted(component))
    triangles=empty_triples=0
    for i in range(len(laws)):
        for j in sorted(w for w in adjacency[i] if w>i):
            for k in sorted(w for w in adjacency[i]&adjacency[j] if w>j):
                triangles+=1;empty_triples+=not bool(pairvalues[i,j]&pairvalues[i,k]&pairvalues[j,k])
    need((len(components),triangles,empty_triples)==(312,26,0),'finite coefficient overlap graph controls')
    out={'scope':'Complete all-positive-parameter equality graph for the540 certified binary coefficient laws. Offsets and dyadic valuations bound every collision search. Counts are of coefficient values, not marked tournament states. Different endpoint sources remain necessary for continuing construction.',
         'laws':laws,'eventual_dyadic_classes':len(classes),'intersecting_law_pairs':len(edges),'collision_values':len(coll),'maximum_collision_value':max(coll) if coll else None,'maximum_value_multiplicity':max(map(len,coll.values())) if coll else 1,
         'isolated_laws':sum(not degrees[i] for i in range(len(laws))),'degree_histogram':dict(Counter(degrees.get(i,0) for i in range(len(laws)))),
         'connected_components':len(components),'component_size_histogram':dict(Counter(map(len,components))),'cycle_rank':len(edges)-len(laws)+len(components),'triangles':triangles,'pairwise_intersecting_triples_without_common_value':empty_triples,
         'asymptotic':'number of distinct supplied coefficients <=X is540*log2(X)+O(1)',
         'exact_dyadic_height_tail':'540*b-6078 for every b>=21','dyadic_tail_from_exponent':threshold,
         'finite_control_height':bound,'finite_control_distinct_values':len(literal),'levels':levels,'edges':edges,'collisions':[{'c':c,'laws':sorted(ids)} for c,ids in sorted(coll.items())]}
    save('coefficient_intersections.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['laws','levels','edges','collisions']}),flush=True)
if __name__=='__main__':run()
