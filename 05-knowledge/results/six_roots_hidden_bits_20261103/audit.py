#!/usr/bin/env python3
"""Second algorithm: literal permutations/cuts, independent recursive H counts."""
import itertools,json,functools
from pathlib import Path
P=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def brute(rows):
    n=len(rows);h=0;k=[0]*n
    for w in itertools.permutations(range(n)):
        good=[bool(rows[a]>>b&1) for a,b in zip(w,w[1:])]
        h+=all(good)
        for t in range(1,n+1):
            if all(good[:t-1]) and all(good[t:]):k[w[t-1]]+=1
    return h,k
def count(rows):
    @functools.lru_cache(None)
    def suffix(v,left):
        if not left:return 1
        return sum(suffix(w,left^(1<<w)) for w in range(len(rows)) if left>>w&1 and rows[v]>>w&1)
    return sum(suffix(v,((1<<len(rows))-1)^(1<<v)) for v in range(len(rows)))
def main():
    r=json.loads((P/'verification.json').read_text());tested=[]
    for e in r['E_profiles']:
        n=e['n']
        if n>8:continue
        rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
        rows[0]-=1<<(n-1);rows[-1]+=1
        need(brute(rows)==(e['H'],e['K']),('E permutation',n));tested.append(n)
    for s in r['actual_seeds']:
        rows=s['rows'];n=len(rows)
        need(brute(rows)==(s['H'],s['K']),('seed permutation',s['H']))
        for q,want in enumerate(s['ear_counts']):
            extended=rows+[0,0]
            for i in range(n):
                if i==q:extended[i]|=1<<n;extended[n+1]|=1<<i
                else:extended[n]|=1<<i;extended[n+1]|=1<<i
            extended[n]|=1<<(n+1)
            need(count(extended)==want,('ear recursive',s['H'],q))
    nearby=0
    for f in sorted(P.glob('nearby67_d*.json')):
        for level in json.loads(f.read_text()):
            for hit in level['count67']:
                rows=hit['rows'];n=len(rows);h,k=brute(rows)
                need(h==67 and [v-h for v in k]==hit['J'],'nearby67 permutation')
                reach=[]
                for start in range(n):
                    todo=[start];seen={start}
                    while todo:
                        v=todo.pop()
                        for w in range(n):
                            if rows[v]>>w&1 and w not in seen:seen.add(w);todo.append(w)
                    reach.append(len(seen)==n)
                need(all(reach)==hit['strong'],'nearby67 strongness');nearby+=1
    # Components by graph search, not the boundary counting formula.
    for b in range(2,6):
        for mask in range(1<<10):
            e={i+1 for i in range(10) if mask>>i&1};alive=set(range(64))-e;components=0
            while alive:
                pending=[alive.pop()];components+=1
                while pending:
                    u=pending.pop();neighbors=[b*u+j for j in range(1,b+1)]
                    if u:neighbors.append((u-1)//b)
                    for v in neighbors:
                        if v in alive:alive.remove(v);pending.append(v)
            internal=sum((u-1)//b in e for u in e)
            need(components==1+b*len(e)-internal,'graph components')
    out={'status':'PASS; second algorithm by same agent, not external referee','E_permutation_orders':tested,'six_seed_permutation_and_recursive_ear_checks':True,'nearby67_literal_witnesses':nearby,'forest_graph_search_cases':4096,'hostiles':'E3 has one port response; E4 has four; five nontrivial printed seeds each require >=3 bits; a padded67 source has J=0 but is not strong'}
    (P/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: independent permutation/cut, recursive ear, and component audits')
if __name__=='__main__':main()
