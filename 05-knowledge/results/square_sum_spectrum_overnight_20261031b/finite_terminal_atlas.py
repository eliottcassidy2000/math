#!/usr/bin/env python3
"""Exact finite affine chart census over binary source words, not a tail proof."""
from collections import defaultdict,deque
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import sys
import sympy as sp

PREVIOUS=Path(__file__).resolve().parent.parent/'square_sum_spectrum_overnight_20261030'
sys.path.insert(0,str(PREVIOUS))
from audit_feedback import gaps,require,subset_dp,tournament


def main():
    # Four is the observable dimension on this source alphabet, not the raw
    # reachable dimension in the unreduced 36-state gap machine.
    representation=json.loads((PREVIOUS/'integral_feedback.json').read_text())
    M=list(map(sp.Matrix,representation['integer_matrices']))
    u=sp.Matrix(representation['initial'])
    prefixes=[(),(1,),(1,0),(1,1)]
    rows=[]
    for p in prefixes:
        row=u
        for a in p:row=row*M[a]
        rows.append(row)
    R=sp.Matrix.vstack(*rows)
    require(R.rank()==4,'source rank')
    require(all(R.col_join(R*M[a]).rank()==4 for a in (0,1)),'source row closure')
    affine=R.row_join(sp.ones(4,1))
    require(affine.rank()==4,'constant coordinate belongs to source representation')
    require(all(affine.col_join(affine*sp.diag(M[a],1)).rank()==4 for a in (0,1)),
            'affine source closure including the constant term')
    terminals=[w for n in range(4) for w in product(range(4),repeat=n)]
    blocks=[w for n in range(3) for w in product((0,1),repeat=n)]
    cache={}
    def column(A):
        if A not in cache:cache[A]=tuple(gaps(p+A) for p in prefixes)
        return cache[A]
    columns=[column(A) for A in terminals]
    edges=[];adj=defaultdict(list);rev=defaultdict(list)
    for i,source in enumerate(columns):
        r=next(j for j in range(1,4) if source[j]!=source[0])
        for B in blocks:
            for j,A in enumerate(terminals):
                target=column(B+A)
                a=Fraction(target[r]-target[0],source[r]-source[0]);b=target[0]-a*source[0]
                if all(target[k]==a*source[k]+b for k in range(4)):
                    require(a>0,'a genuine universal map has positive slope')
                    index=len(edges);edges.append((i,j,B,a,b));adj[i].append(index);rev[j].append(index)
    # Reachability gives an independently transparent SCC partition for 85 nodes.
    def reachable(i):
        seen={i};q=[i]
        for x in q:
            for e in adj[x]:
                y=edges[e][1]
                if y not in seen:seen.add(y);q.append(y)
        return seen
    reach=[reachable(i) for i in range(len(terminals))]
    components=[];owner={}
    for i in range(len(terminals)):
        if i in owner:continue
        c=sorted(j for j in reach[i] if i in reach[j])
        for j in c:owner[j]=len(components)
        components.append(c)
    # In each SCC assign affine gauges from one root; every edge must have the
    # identical transition between gauges. This checks cycle identity without
    # enumerating exponentially many directed cycles.
    cyclic_edges=0;nonempty_cyclic_edges=0
    for component in components:
        root=component[0];gauge={root:(Fraction(1),Fraction(0))};q=[root]
        for i in q:
            for e in adj[i]:
                _,j,B,a,b=edges[e]
                if owner[i]!=owner[j]:continue
                current=(a*gauge[i][0],a*gauge[i][1]+b)
                if j not in gauge:gauge[j]=current;q.append(j)
                require(gauge[j]==current,'nontrivial affine cycle')
                cyclic_edges+=1
                if B:
                    nonempty_cyclic_edges+=1
                    require(set(B)<={1},'productive zero on a closed chart return')
        require(len(gauge)==len(component),'SCC gauge coverage')
    # Actual finite reachability in (chart, count), retaining all distinct counts.
    # It terminates here; the all-size theorem separately proves finiteness.
    states={(i,gaps(A)) for i,A in enumerate(terminals)};q=deque(states)
    while q:
        i,h=q.popleft()
        for e in adj[i]:
            _,j,B,a,b=edges[e];k=a*h+b
            require(k.denominator==1 and k>0,'actual chart count integrality')
            state=(j,int(k))
            if state not in states:states.add(state);q.append(state)
        require(len(states)<1000000,'finite exploration cap')
    # Independent actual tournaments for every discovered relation at four
    # short source words, beyond the matrix-only closure proof.
    controls=0
    for i,j,B,a,b in edges:
        for w in ((),(0,),(1,),(1,0)):
            require(subset_dp(tournament(w+B+terminals[j]))==a*subset_dp(tournament(w+terminals[i]))+b,
                    'literal affine edge')
            controls+=1
    examples=[]
    for i,j,B,a,b in edges:
        if a!=1 and B and len(examples)<10:
            examples.append({'source_terminal':terminals[i],'inserted':B,'target_terminal':terminals[j],
                             'a':str(a),'b':str(b)})
    singleton_controls=0
    for m in range(15):
        w=(1,)*m
        require(gaps(w+(1,1,0))==2*gaps(w+(1,0))-1,'one-letter doubling hostile')
        if m<=4:
            require(subset_dp(tournament(w+(1,1,0)))==2*subset_dp(tournament(w+(1,0)))-1,
                    'literal singleton hostile')
        singleton_controls+=1
    print(json.dumps({'status':'passed','source_alphabet':[0,1],'source_dimension':4,
                      'affine_source_dimension':4,'singleton_hostile_controls':singleton_controls,
                      'prefix_basis':prefixes,'terminal_length_at_most':3,'block_length_at_most':2,
                      'terminal_charts':len(terminals),'affine_edges':len(edges),
                      'scc_sizes':sorted(map(len,components),reverse=True),'within_scc_edges':cyclic_edges,
                      'nonempty_within_scc_edges':nonempty_cyclic_edges,
                      'reachable_chart_count_pairs':len(states),'reachable_H_values':sorted({h for _,h in states}),
                      'independent_subset_dp_controls':controls,'nonunit_slope_examples':examples,
                      'scope':'Complete only for this finite atlas; all-size no-pumping follows separately from the fixed-terminal rigidity theorem.'},indent=2,sort_keys=True))


if __name__=='__main__':main()
