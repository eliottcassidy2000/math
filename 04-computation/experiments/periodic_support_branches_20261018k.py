#!/usr/bin/env python3
"""Decode the complete 375-branch necessary cover of a six-way contact.

Keep original factor vertices. Congruences are difference equations modulo 2
or 3, not colour equalities or a quotient of the geometric source.
"""
import argparse
import collections
import json
from pathlib import Path


def cosets(D):
    return [(p,next(iter({d%p for d in D}))) for p in (2,3) if len({d%p for d in D})==1]


def solve_voltages(n,edges,p):
    """Return vertex potentials including isolates, or an oriented bad cycle."""
    assert p in (2,3)
    adj=[[] for _ in range(n)];tree=[[] for _ in range(n)]
    for u,v,r in edges:
        assert 0<=u<n and 0<=v<n
        adj[u].append((v,r%p));adj[v].append((u,-r%p))
    potential=[None]*n;components=0
    for root in range(n):
        if potential[root] is not None:continue
        components+=1;potential[root]=0;queue=collections.deque([root])
        while queue:
            u=queue.popleft()
            for v,r in adj[u]:
                expect=(potential[u]+r)%p
                if potential[v] is None:
                    potential[v]=expect;queue.append(v)
                    tree[u].append((v,r));tree[v].append((u,-r%p))
                elif potential[v]!=expect:
                    # The edge u->v plus the already assigned tree path v->u.
                    parent={v:None};q=collections.deque([v])
                    while u not in parent:
                        a=q.popleft()
                        for b,s in tree[a]:
                            if b not in parent:parent[b]=(a,s);q.append(b)
                    path=[];a=u
                    while a!=v:
                        b,s=parent[a];path.append((b,a,s));a=b
                    cycle=[(u,v,r)]+list(reversed(path))
                    assert all(cycle[i][1]==cycle[(i+1)%len(cycle)][0] for i in range(len(cycle)))
                    assert sum(e[2] for e in cycle)%p!=0
                    allowed={(a,b,t%p) for a,b,t in edges}|{(b,a,-t%p) for a,b,t in edges}
                    assert all(e in allowed for e in cycle)
                    return dict(consistent=False,modulus=p,cycle=cycle,cycle_sum=sum(e[2] for e in cycle)%p)
    assert all((potential[v]-potential[u]-r)%p==0 for u,v,r in edges)
    return dict(consistent=True,modulus=p,components=components,potential=potential)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    masks=0
    for m in range(1,64):
        D={d for d in range(6) if m>>d&1};rich=any((d+1)%6 in D for d in D)
        assert bool(cosets(D))==(not rich);masks+=1
    all_branches={(i,p,a,j,q,b) for i in range(6) for j in range(i+1,6)
                  for p in (2,3) for a in range(p) for q in (2,3) for b in range(q)}
    assert len(all_branches)==375
    controls=[solve_voltages(4,[(0,1,1),(1,2,1),(2,0,1)],2),
              solve_voltages(4,[(0,1,1),(1,2,1),(2,0,1)],3),
              solve_voltages(4,[(0,1,1),(1,2,1),(2,0,0)],3)]
    assert [c['consistent'] for c in controls]==[False,True,False]
    assert controls[1]['components']==2 # the isolated original vertex survives
    data=json.loads(args.input.read_text())
    groups=[g for g in data['support_groups'] if g['family']=='u_times_v_square']
    common=[None]*len(groups);frequencies=[collections.Counter() for _ in groups]
    rows=[];linear_systems=equations=0
    for w in data['witnesses']:
        selected=[];systems={(axis,p):[] for axis in range(6) for p in (2,3)}
        for gi,g in enumerate(groups):
            options=[]
            for axis,ps in g['factors']:
                D={(w['colours'][axis][j]-w['colours'][axis][i])%6 for i,j in ps}
                options.append(cosets(D))
            branches={(i,p,a,j,q,b) for i in range(6) for j in range(i+1,6)
                      for p,a in options[i] for q,b in options[j]}
            assert branches and branches<=all_branches
            common[gi]=branches if common[gi] is None else common[gi]&branches
            frequencies[gi].update(branches)
            branch=min(branches);selected.append(branch)
            i,p,a,j,q,b=branch
            for pos,mod,res in ((i,p,a),(j,q,b)):
                axis,ps=g['factors'][pos]
                systems[axis,mod].extend((x,y,res) for x,y in ps)
        components=[]
        for (axis,p),edges in sorted(systems.items()):
            answer=solve_voltages(data['coordinate_sizes'][axis],edges,p)
            assert answer['consistent'];linear_systems+=1;equations+=len(edges)
            assert all((w['colours'][axis][j]-w['colours'][axis][i]-r)%p==0 for i,j,r in edges)
            components.append([axis,p,answer['components']])
        rows.append(dict(pattern=w['pattern'],selected_branches=selected,component_counts=components))
    summaries=[]
    for gi,(inter,counts) in enumerate(zip(common,frequencies)):
        top,count=counts.most_common(1)[0]
        summaries.append(dict(group=gi,common_branches=sorted(inter),most_frequent_branch=top,
                              most_frequent_count=count,observed_branches=len(counts)))
    out=dict(status='FINITE-EXACT; branch-cover theorem is independent of these witnesses',
             support_masks_checked=masks,branches_per_contact=len(all_branches),controls=controls,
             patterns=len(rows),six_way_groups=len(groups),linear_systems=linear_systems,equations=equations,
             group_summaries=summaries,source_preserving_decodings=rows,
             scope='Necessary cover for additive colourings; congruence consistency alone is not a graph-colouring certificate.')
    args.out.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print({k:v for k,v in out.items() if k!='source_preserving_decodings'})


if __name__=='__main__':main()
