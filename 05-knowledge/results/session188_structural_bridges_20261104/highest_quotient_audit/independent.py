#!/usr/bin/env python3
"""Independent highest-change audit via original-module edge-code spans."""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'covered_laminar'))
from verify_ordered import tournament,induced,modules,profile,hp,extend,need
ROOT=Path(__file__).resolve().parent.parent


def ordered_cover_profiles(T,U,children):
    return all(profile(induced(T,[i for i in range(len(T)) if A>>i&1]))==profile(induced(U,[i for i in range(len(U)) if A>>i&1])) for A in children)


def original_tree(T):
    n=len(T);mods=modules(T);allmods=mods+[1<<i for i in range(n)]
    strong=[A for A in allmods if not any(A&B and A&B!=A and A&B!=B for B in allmods)]
    result={}
    for A in strong:
        if A.bit_count()<2:continue
        proper=[B for B in strong if B!=A and A&B==B]
        children=[B for B in proper if not any(B!=C and B&C==B for C in proper)]
        need(sum(children)==A,'children partition')
        result[A]=children
    return mods,result


def Q(T,children):
    reps=[(a&-a).bit_length()-1 for a in children]
    return induced(T,reps)


def is_linear(Q):return len(set(r.bit_count() for r in Q))==len(Q)


def tops_for(T,U,tree):
    changed=[A for A,C in tree.items() if Q(T,C)!=Q(U,C)]
    tops=[A for A in changed if not any(A!=B and A&B==A for B in changed)]
    return tops,changed


def main():
    totals=dict(families=0,equal_families=0,equal_nonidentity=0,linear_highest_nodes=0,prime_highest_nodes=0,multiple_highest_families=0,prime_only_families=0)
    queries=0;literal=0;byorder=[]
    for n in range(1,6):
        saved=json.loads((ROOT/'simultaneous_response'/f'full_response_n{n}.json').read_text());basis=saved['basis'];records=saved['records'];graphs=[tournament(n,k) for k in range(len(records))];vectors=[r['vector'] for r in records];H=[v[0] for v in vectors];pairs=list(combinations(range(n),2));local={k:0 for k in totals}
        for source,T in enumerate(graphs):
            mods,tree=original_tree(T);span={0}
            for A in mods:
                x=sum(1<<j for j,(u,v) in enumerate(pairs) if A>>u&1 and A>>v&1)
                nxt={s^x for s in span};need(not span&nxt,'module incidence independent');span|=nxt
            for x in span:
                target=source^x;local['families']+=1
                if H[source]!=H[target]:continue
                local['equal_families']+=1;U=graphs[target];tops,changed=tops_for(T,U,tree)
                need(bool(tops)==bool(x),'highest existence')
                need(all(not A&B for A,B in combinations(tops,2)),'highest antichain')
                for bit,(u,v) in enumerate(pairs):
                    if x>>bit&1:need(any(A>>u&1 and A>>v&1 for A in tops),'changed edge coverage')
                if not tops:continue
                local['equal_nonidentity']+=1;lin=0
                for A in tops:
                    children=tree[A];q,r=Q(T,children),Q(U,children)
                    need(ordered_cover_profiles(T,U,children),'complete child profiles equal')
                    if not is_linear(q):
                        need(not any(m!=(1<<len(q))-1 for m in modules(q)),'highest non-linear quotient prime')
                        need(all(r[i]==(((1<<len(q))-1)^q[i]^(1<<i)) for i in range(len(q))),'prime converse')
                        local['prime_highest_nodes']+=1;continue
                    need(is_linear(r),'linear remains linear');lin+=1;local['linear_highest_nodes']+=1
                    old=sorted(range(len(q)),key=lambda i:-q[i].bit_count());new=sorted(range(len(r)),key=lambda i:-r[i].bit_count())
                    cut=next(t for t in range(1,len(q)) if set(old[:t])!=set(new[:t]));prefix=sum(children[i] for i in old[:cut]);rep=(A&-A).bit_length()-1
                    outside=sum(1<<u for u in range(n) if not A>>u&1 and T[u]>>rep&1);a,c=outside,outside|prefix
                    def value(v):return sum(z for z,(s,t) in zip(v,basis) if s&a==s and t&c==t)
                    f,g=value(vectors[source]),value(vectors[target]);ambient=Q(T,[A]+[1<<u for u in range(n) if not A>>u&1]);bound=2*hp(ambient)
                    for B in children:bound*=hp(induced(T,[u for u in range(n) if B>>u&1]))
                    need(g-f>=bound,'literal physical cut lower bound');queries+=2
                    # A separately computed literal integer DP for every distinct source/target/cut.
                    need((f,g)==(hp(extend(T,a,c)),hp(extend(U,a,c))),'independent full response DP');literal+=2
                local['multiple_highest_families']+=len(tops)>1;local['prime_only_families']+=lin==0
        for k,v in local.items():totals[k]+=v
        byorder.append(dict(n=n,**local))
    produced=json.loads((ROOT/'highest_quotient_reduction'/'normal.json').read_text())
    need(totals==produced['totals'] and byorder==produced['by_order'],'complete producer census')
    controls=[]
    for name in ['multiple_highest_control','prime_three_control']:
        d=produced[name];T,U=tuple(d['T']),tuple(d['U']);_,tree=original_tree(T);tops,changed=tops_for(T,U,tree)
        if name=='multiple_highest_control':need(len(tops)==2 and sorted(len(tree[A]) for A in tops)==[2,3],'linear plus prime antichain')
        else:need(tops==[63] and len(changed)==4,'prime root shields three changed linear descendants')
        controls.append(dict(name=name,highest=tops,all_changed=changed))
    return dict(status='PASS',algorithm='Original-module incidence spans, explicit original strong trees, complete response-bank queries and independent literal integer DP.',totals=totals,by_order=byorder,full_response_query_comparisons=queries,literal_integer_dp_comparisons=literal,controls=controls)
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
