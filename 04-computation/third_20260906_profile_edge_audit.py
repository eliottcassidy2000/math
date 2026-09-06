#!/usr/bin/env python3
"""Independent scalar-projection-first audit of the selected sheet edge.

No repository producer imports. Full words are consulted only after the
complete connected scalar universe is reconstructed. The actual control
uses modulo-t and amplitude proofs for all mixed relations.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import gcd, prod
from pathlib import Path
import json
import sys
import sympy as S

ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
Q=91**6
GATES=0
DIGEST=sha256()

def gate(ok,label,data=None):
    global GATES
    if not ok:raise RuntimeError(f'{label}: {data!r}')
    GATES+=1
    DIGEST.update((label+':'+repr(data)+'\n').encode())

def load():
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    gate(sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','pinned-inherited-profile-bytes')
    levels=json.loads(raw)['levels']
    profiles={int(k):{(c,tuple(w)) for c,w in entry['profiles']} for k,entry in levels.items()}
    allowed={k:{c for c,w in rows} for k,rows in profiles.items()}
    for k,entry in levels.items():
        gate(allowed[int(k)]==set(entry['gcds']),'exact-scalar-projection',int(k))
        gate(all(len(w)==int(k) and tuple(sorted(w))==w for c,w in profiles[int(k)]),'profile-types',int(k))
    return profiles,allowed

def components(n,edges):
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in edges:
        a,b=find(i),find(j)
        if a!=b:parent[b]=a
    return sorted(tuple(i for i in range(n) if find(i)==r) for r in {find(i) for i in range(n)})

def full_failure(D,profiles):
    for size in range(1,7):
        for I in combinations(range(7),size):
            c=gcd(*(D[i] for i in I))
            word=tuple(sorted(gcd(c,D[j]) for j in range(7) if j not in I))
            if (c,word) not in profiles[7-size]:return I,c,word
    return None

def enumerate_quotient(profiles,allowed,threshold):
    values=sorted(v for v in allowed[6] if v>=threshold)
    adjacency={v:{w for w in values if gcd(v,w)>=threshold} for v in values}
    # The producer prunes by partial FULL WORDS. This audit uses only scalar
    # projections until the complete seven-state bank has been generated.
    @lru_cache(None)
    def scalar_valid(D):
        for size in range(2,len(D)+1):
            possible=allowed[7-size] if size<=6 else {1}
            if any(gcd(*I) not in possible for I in combinations(D,size)):return False
        return True
    current={(v,) for v in values};counts=[len(current)]
    for size in range(2,8):
        extensions=set()
        for D in sorted(current):
            for v in sorted(set().union(*(adjacency[d] for d in D))):
                candidate=tuple(sorted(D+(v,)))
                if scalar_valid(candidate):extensions.add(candidate)
        current=extensions;counts.append(len(current))
    expected={7:[36,104,377,1031,1923,1700,104],6:[37,144,788,2929,6407,5985,493]}
    gate(counts==expected[threshold],'independent-scalar-only-layer-counts',(threshold,counts))
    survivors=[];rejected=Counter()
    for D in sorted(current):
        edges=[(i,j) for i,j in combinations(range(7),2) if gcd(D[i],D[j])>=threshold]
        gate(components(7,edges)==[tuple(range(7))],'complete-state-connectivity',(threshold,D))
        failure=full_failure(D,profiles)
        if failure is None:
            survivors.append(D)
            gate(gcd(*D)==1,'complete-full-profile-survivor',D)
        else:
            I,c,word=failure
            gate((c,word) not in profiles[7-len(I)],'literal-final-word-rejection',(threshold,D,failure))
            rejected[len(I)]+=1
    gate(len(survivors)==(0 if threshold==7 else 26),'complete-profile-survivor-count',threshold)
    if threshold==7:
        hostile=(8,8,9,9,10,60,72)
        gate(hostile in current and full_failure(hostile,profiles)==((5,),60,(3,3,4,4,10,12)),
             'exact-scalar-projection-hostile')
        print('EXACT_SCALAR_HOSTILE',hostile,'excluded',(60,(3,3,4,4,10,12)))
    print('INDEPENDENT_QUOTIENT',json.dumps({'threshold':threshold,'scalar_layer_counts':counts,
        'rejection_subset_sizes':dict(sorted(rejected.items())),'survivors':survivors},sort_keys=True,separators=(',',':')))
    return set(survivors)

def strict_atlas():
    # Generate allowed primitive sums multiplicatively from a prime sieve,
    # independently of the producer's per-integer factorization predicate.
    sieve=[True]*357;sieve[0]=sieve[1]=False
    for p in range(2,19):
        if sieve[p]:
            for m in range(p*p,357,p):sieve[m]=False
    primes=[p for p in range(2,357) if sieve[p] and p%3==2]
    sums={1}
    for p in primes:
        sums={v*p**e for v in sums for e in range(3) if v*p**e<=356}
    sums.discard(1)
    gate(25 in sums and 125 not in sums,'strict-cube-exponent-boundary')
    atlas={(a,s-a) for s in sums for a in range(1,(s+1)//2) if gcd(a,s-a)==1}
    gate(len(atlas)==5855,'independently-generated-strict-atlas')
    return atlas

def actual_control(profiles,survivors):
    V=(1,2,3,4,5,6);U=(12584,14872,117,9999,98890,132990,10296)
    t=360*(1000*Q+1);row=tuple(t*v for v in V)+U
    states=tuple(gcd(t,u) for u in U)
    gate(states==(8,8,9,9,10,30,72) and states in survivors,'actual-sheet-realization')
    gate(gcd(*V)==gcd(*U)==gcd(*row)==1 and len(set(row))==13 and min(row)>0,'actual-positive-distinct-primitive')
    gate(gcd(1000*Q+1,prod(U))==1,'actual-scale-sidecar')
    gate(t==204432930734760360 and sum(row)==4293091545430247308 and sum(row)<Q*Q,'actual-physical-box')
    atlas=strict_atlas();edges=[];relations=[]
    for i,j in combinations(range(13),2):
        common=gcd(row[i],row[j]);a,b=row[i]//common,row[j]//common
        if tuple(sorted((a,b))) in atlas:
            edges.append((i,j));r=[0]*13;r[i]=b;r[j]=-a;relations.append(r)
            gate(sum(v*w for v,w in zip(r,row))==0 and max(map(abs,r))<=355,'literal-edge-row',(i,j))
    gate(components(13,edges)==[tuple(range(6)),tuple(range(6,13))],'actual-exact-decoder-components')
    gate(S.Matrix(relations).rank()==11,'independent-exact-rank11')
    tail_edges=[]
    for i,j in edges:
        if i>=6:
            a,b=i-6,j-6;d=gcd(U[a],U[b])
            tail_edges.append((a,b,U[a]//d,U[b]//d,gcd(t,d)))
    expected=[(0,6,11,9,8),(1,6,13,9,8),(2,6,1,88,9),(3,6,101,104,9),(4,5,29,39,10),(5,6,155,12,6)]
    gate(tail_edges==expected and min(e[-1] for e in tail_edges)==6,'actual-edge-bound-attainment')
    # All mixed signed relations, using independent global inequalities.
    # t | c*u forces t/gcd(t,u) | c in the first orientation.
    mixed=0
    for i,j in combinations(range(6),2):
        for k in range(7):
            gate(t//gcd(t,U[k])>Q,'modulo-t-excludes-mixed-relation',(i,j,k));mixed+=1
    for i,j in combinations(range(7),2):
        for k in range(6):
            gate(t*V[k]>Q*(U[i]+U[j]),'amplitude-excludes-reverse-mixed-relation',(i,j,k));mixed+=1
    gate(mixed==231,'complete-mixed-support-universe')
    for component in (V,U):
        for a,b in combinations(component,2):
            gate(max(a,b)//gcd(a,b)<=Q,'actual-internal-pair-height',(a,b))
    profile_count=0;maxima=[]
    for size in range(7,13):
        maximum=0
        for I in combinations(range(13),size):
            c=gcd(*(row[i] for i in I))
            w=tuple(sorted(gcd(c,row[j]) for j in range(13) if j not in I))
            gate((c,w) in profiles[13-size],'all-physical-profile-words',(I,c,w))
            maximum=max(maximum,c);profile_count+=1
        maxima.append(maximum)
    gate(profile_count==4095 and maxima==[72,10,9,3,2,1],'complete-physical-profile-universe')
    clearance=Fraction(min(min(v%7,(-v)%7) for v in row),7)
    gate(clearance==Fraction(1,7),'actual-safe-phase')
    print('INDEPENDENT_ACTUAL_CONTROL',json.dumps({'V':V,'U':U,'t':t,'sheet_states':states,
        'sum':sum(row),'edges':edges,'tail_edges':tail_edges,'rank':11,'mixed_supports':mixed,
        'profiles':profile_count,'max_subset_gcds':maxima,'clearance':str(clearance)},sort_keys=True,separators=(',',':')))

def main():
    profiles,allowed=load()
    enumerate_quotient(profiles,allowed,7)
    survivors=enumerate_quotient(profiles,allowed,6)
    producer=(ROOT/'05-knowledge/results/third_20260906_decoder_profile_graph.out').read_text()
    line=next(line for line in producer.splitlines() if line.startswith('EXHAUSTIVE_PROFILE_GRAPH threshold=6 '))
    producer_survivors={tuple(row) for row in json.loads(line.split('survivors=',1)[1])}
    gate(survivors==producer_survivors,'all26-independent-producer-survivor-comparison')
    actual_control(profiles,survivors)
    print('PASS',GATES,'always-active exact gates')
    print('SEMANTIC_SHA256',DIGEST.hexdigest())

if __name__=='__main__':main()
