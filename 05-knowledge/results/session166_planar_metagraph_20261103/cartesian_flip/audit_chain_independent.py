"""Independent root/diamond enumeration, separate from generic induced matcher."""
import json
from itertools import combinations
from pathlib import Path
BASE=Path(__file__).resolve().parent

def search(order,m):
    host=json.loads((BASE/f'metagraph_n{order}.json').read_text())
    n=len(host['representatives']);adj=[set()for _ in range(n)]
    for a,b in host['edges']:adj[a].add(b);adj[b].add(a)
    calls=0;diamonds=0
    for root in range(n):
        candidates=[]
        for a,b in combinations(sorted(adj[root]),2):
            if b not in adj[a]:continue
            for t in sorted((adj[a]&adj[b])-adj[root]-{root}):
                candidates.append((t,a,b))
        diamonds+=len(candidates)
        def dfs(chosen,used):
            nonlocal calls
            calls+=1
            if len(chosen)==m:return chosen
            for c in candidates:
                if any(x in used for x in c):continue
                if chosen and c[0] not in adj[chosen[-1][0]]:continue
                valid=True
                for k,old in enumerate(chosen):
                    for i,x in enumerate(old):
                        for j,y in enumerate(c):
                            expected=(k==len(chosen)-1 and i==j==0)
                            if (y in adj[x])!=expected:valid=False;break
                        if not valid:break
                    if not valid:break
                if valid:
                    got=dfs(chosen+[c],used|set(c))
                    if got is not None:return got
            return None
        got=dfs([],{root})
        if got is not None:return {'order':order,'m':m,'found':[root]+[x for c in got for x in c],'calls':calls,'diamonds_seen':diamonds}
    return {'order':order,'m':m,'found':None,'calls':calls,'diamonds_seen':diamonds}

if __name__=='__main__':
    out=[]
    for n,m,expected in [(5,2,False),(6,3,False),(7,4,False),(6,2,True),(7,3,True)]:
        item=search(n,m)
        if (item['found']is not None)!=expected:raise RuntimeError('independent chain search regression')
        out.append(item);print(json.dumps(item),flush=True)
    (BASE/'audit_chain_independent.json').write_text(json.dumps(out,indent=2)+'\n')
