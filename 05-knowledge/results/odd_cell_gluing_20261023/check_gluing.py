"""Exhaustive tiny controls for partial-bijection sections and chart-conditioned flow."""
from itertools import combinations,permutations,product
from collections import Counter
from pathlib import Path
import argparse,json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)

def sections(domains,relations):
    """Connected base, actual partial bijections, including empty relations."""
    n=len(domains);V={(i,x) for i,D in enumerate(domains) for x in D}
    adj={v:set() for v in V};deg=[0]*n
    for (i,j),R in relations.items():
        if len(set(R.values()))!=len(R):raise ValueError('not a partial bijection')
        deg[i]+=1;deg[j]+=1
        for x,y in R.items():
            if (i,x) not in V or (j,y) not in V:raise ValueError('state outside declared domain')
            adj[i,x].add((j,y));adj[j,y].add((i,x))
    seen=set();out=[]
    for v in sorted(V):
        if v in seen:continue
        seen.add(v);C=[v]
        for w in C:
            for z in adj[w]:
                if z not in seen:seen.add(z);C.append(z)
        # Completeness matters: coverage alone can miss an undefined chord.
        if len(C)==n and len({i for i,x in C})==n and all(len(adj[i,x])==deg[i] for i,x in C):
            out.append(tuple(x for i,x in sorted(C)))
    return sorted(out)

def main():
    partial=[]
    for k in range(4):
        for U in combinations(range(3),k):
            for V in permutations(range(3),k):partial.append(dict(zip(U,V)))
    need(len(partial)==34,'partial bijection count')
    tally=Counter()
    for r,s,t in product(partial,repeat=3):
        relations={(0,1):r,(1,2):s,(0,2):t}
        truth=[(a,b,c) for a,b,c in product(range(3),repeat=3) if r.get(a)==b and s.get(b)==c and t.get(a)==c]
        got=sections([range(3)]*3,relations);need(got==truth,'state graph changed a joint choice')
        tally[len(got)]+=1
    # Base has connected state support spanning its cells, but the chord is undefined.
    hole={(0,1):{0:0},(1,2):{0:0},(0,2):{}}
    need(sections([[0]]*3,hole)==[],'forgot missing chord domain')
    # Private path suppression composes actual partial maps exactly.
    eliminated=0
    for r,s in product(partial,repeat=2):
        composite={x:s[y] for x,y in r.items() if y in s}
        direct={(x,z) for x,z in product(range(3),repeat=2) if any(r.get(x)==y and s.get(y)==z for y in range(3))}
        need(set(composite.items())==direct,'private path composition');eliminated+=1
    # Correlated chart branches: each block may swap only with the SAME bit.
    # P counts first-block off-diagonal entries minus second-block ones.
    alternatives=[]
    for a,b in product(range(2),repeat=2):
        X=[[0]*4 for _ in range(4)]
        for k,t in enumerate((a,b)):
            for i in range(2):X[2*k+i][2*k+(i^t)]=1
        P=sum(X[i][j] for i,j in [(0,1),(1,0)])-sum(X[i][j] for i,j in [(2,3),(3,2)])
        alternatives.append({'first_bit':a,'second_bit':b,'projected_cost':P,'shared_chart':a==b})
    need({x['projected_cost'] for x in alternatives if x['shared_chart']}=={0},'joint branches disagree')
    need({x['projected_cost'] for x in alternatives}=={-2,0,2},'union should add phantom projection')
    return {'status':'FINITE-EXACT','triangle_partial_bijection_systems':34**3,
            'counts_by_number_of_sections':sorted(tally.items()),'private_path_compositions':eliminated,
            'conditioned_transport_hostile':alternatives,'checks':CHECKS}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=json.loads(json.dumps(main()));f=HERE/'gluing.json'
    if a.write:f.write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads(f.read_text())==out,'receipt mismatch')
    print(json.dumps(out,indent=2))
