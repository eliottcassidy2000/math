#!/usr/bin/env python3
"""Independent graph/component and matrix audit; imports no producer code."""
from collections import defaultdict
from itertools import permutations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def require(ok,why):
    if not ok:raise RuntimeError(why)

def components_pairing(blocks,flip):
    # Inspect entire undirected multigraph, not the producer's alternating walk.
    owner={x:c for c,xs in blocks.items() for x in xs}
    representatives={min(xs):c for c,xs in blocks.items()}
    require(len(owner)==sum(map(len,blocks.values())),'overlapping blocks')
    adjacency={x:[flip(x)] for x in owner}
    require(all(flip(x) in owner and flip(x)!=x and flip(flip(x))==x for x in owner),'tau not free')
    for c,xs in blocks.items():
        require(len(xs)%2==1,'even block')
        seq=sorted(xs)[1:]
        for i in range(0,len(seq),2):
            a,b=seq[i:i+2];adjacency[a].append(b);adjacency[b].append(a)
    seen=set();pairs=set();cycle_components=0
    for start in owner:
        if start in seen:continue
        todo=[start];component=set()
        while todo:
            x=todo.pop()
            if x in component:continue
            component.add(x);todo.extend(adjacency[x])
        seen|=component
        endpoints=[x for x in component if len(adjacency[x])==1]
        require(set(endpoints)==component&set(representatives),'wrong endpoints')
        if not endpoints:
            cycle_components+=1
        else:
            require(len(endpoints)==2,'component not path or cycle')
            pair=tuple(sorted(representatives[x] for x in endpoints))
            require(pair[0]!=pair[1],'same class at both endpoints')
            pairs.add(pair)
    require(len(pairs)*2==len(blocks),'not every class paired')
    return pairs,cycle_components

def odd_partitions_set(m):
    blocks=[]
    def visit(k):
        if k==m:
            if all(len(b)%2 for b in blocks):yield [b[:] for b in blocks]
            return
        for b in blocks:
            b.append(k);yield from visit(k+1);b.pop()
        blocks.append([k]);yield from visit(k+1);blocks.pop()
    yield from visit(0)

def main():
    receipt=json.loads((HERE/'parity_receipt.json').read_text())
    joint=json.loads((HERE/'joint_audit.json').read_text())
    comparisons=[];actual=[]
    for p,j in zip(receipt['census'],joint['census']):
        require(p['n']==j['n'],'order mismatch')
        for a,b in [('classes','ordinary'),('self_converse','self_converse'),('F_n','F'),('SC_H','self_converse_F')]:
            require(p[a]==j[b],f'joint census mismatch {p["n"]} {a}')
        comparisons.append(p['n'])
        n=p['n']
        if n<3:continue
        blocks={int(c):xs for c,xs in p['blue_fibres'].items()}
        q=(n-1)*(n-2)//2;flip=lambda x:x^((1<<q)-1)
        expected,cycles=components_pairing(blocks,flip)
        got={tuple(sorted(item['classes'])) for item in p['refined_pairs'] if item['self_converse']}
        require(got==expected,'component endpoints differ from producer pairing')
        chords=[(i,j) for i in range(n) for j in range(i+2,n)]
        perms=list(permutations(range(n)))
        def graph(mask):
            a=[[int(i<j) for j in range(n)] for i in range(n)]
            for k,(i,j) in enumerate(chords):
                a[i][j]=1-((mask>>k)&1);a[j][i]=(mask>>k)&1
            return a
        def code(a,perm):
            return tuple(a[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        def canon(a):return min(code(a,perm) for perm in perms)
        canonical_classes={}
        hvalues={}
        for c,xs in blocks.items():
            keys={canon(graph(x)) for x in xs}
            require(len(keys)==1,'rho block not one ordinary class')
            key=keys.pop();require(key not in canonical_classes,'two blocks same class')
            canonical_classes[key]=c
            a=graph(min(xs));op=[[a[j][i] for j in range(n)] for i in range(n)]
            require(canon(op)==key,'blue class not self-converse')
            h=sum(all(a[perm[i]][perm[i+1]] for i in range(n-1)) for perm in perms)
            hvalues[c]=h
            for x in xs:
                a=graph(x)
                require(all(a[i][j]==a[n-1-j][n-1-i] for i in range(n) for j in range(n) if i!=j),'not sigma fixed')
        for trail in p['blue_trails']:
            require(trail['H']==[hvalues[c] for c in trail['classes']],'diagram H labels wrong')
        actual.append({'n':n,'self_converse_classes':len(blocks),'blue_masks':sum(map(len,blocks.values())),
                       'class_pairs':len(expected),'discarded_cycle_components':cycles})
    abstract=[]
    for m in range(2,11,2):
        count=0;cycles=0
        for blocks in odd_partitions_set(m):
            _,c=components_pairing(dict(enumerate(blocks)),lambda x:x^1)
            count+=1;cycles+=c
        abstract.append({'presentation_size':m,'all_odd_fibre_partitions':count,'cycle_components':cycles})
    result={'status':'PASS','scope':'independent matrix and component audit; no producer imports',
            'ordinary_joint_census_agrees_n':comparisons,'blue_component_checks':actual,
            'exhaustive_abstract_odd_partitions':abstract}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    (HERE/'audit_parity.json').write_text(text);print(text,end='')

if __name__=='__main__':main()
