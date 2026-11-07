#!/usr/bin/env python3
"""Independent explicit word DFS, with no memoization or pruning except visited vertices."""
import json, math
from pathlib import Path

def paths(n):
    adj={i: [j for j in range(1,n+1) if j!=i and math.isqrt(i+j)**2==i+j] for i in range(1,n+1)}
    result=[]
    def walk(p,used):
        if len(p)==n:
            if p[0]<p[-1]: result.append(p[:])
            return
        for j in adj[p[-1]]:
            if not used & (1<<(j-1)):
                walk(p+[j],used | 1<<(j-1))
    for i in range(1,n+1): walk([i],1<<(i-1))
    return {'n':n,'edge_count':sum(map(len,adj.values()))//2, 'paths_up_to_reversal':result,'count':len(result)}
if __name__=='__main__':
    rows=[paths(n) for n in range(2,26)]
    out={'status':'FINITE-EXACT','universe':'Q_n with vertices 1..n, edges when distinct vertex sum is square; 2 <= n <= 25','method':'Explicit simple-path DFS from every vertex, terminal first<last reversal quotient, no pruning','rows':rows}
    expected={15:1,16:1,17:1,23:3,25:10}
    for r in rows:
        if r['count'] != expected.get(r['n'],0): raise ValueError(r)
    print(json.dumps(out,indent=2))
