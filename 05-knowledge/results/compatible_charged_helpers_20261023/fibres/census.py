"""Complete fixed-margin/full-charge fibres of the authenticated28source menu."""
import json,sys
from collections import defaultdict,Counter
from itertools import combinations
from pathlib import Path
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
def load():
    d=json.loads((RESULTS/'native_packet_decomposition_20261023/native/domains.json').read_text())
    return sorted(d['sources'],key=lambda r:r['source'])
def groups(k,rows):
    X=[r['code'][0] for r in rows];Y=[sum(b<<j for j,b in enumerate(r['code'][1:])) for r in rows]
    out=defaultdict(list)
    for ids in combinations(range(len(rows)),k):
        x=0;y=0;m=0
        for i in ids:x+=X[i];y^=Y[i];m|=1<<i
        out[(x%128,y)].append(m)
    return out
def components(states,budget=6):
    parent=list(range(len(states)))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    active=0
    for i,a in enumerate(states):
        active+=1
        for j in range(i):
            if (a^states[j]).bit_count()<=budget:
                u,v=root(i),root(j)
                if u!=v:parent[u]=v;active-=1
                if active==1:break
    out=defaultdict(list)
    for i,m in enumerate(states):out[root(i)].append(m)
    return sorted((sorted(v) for v in out.values()),key=lambda x:(len(x),x))
def census(k):
    rows=load();g=groups(k,rows);bad=[];hist=Counter();sizes=Counter()
    for key,ss in sorted(g.items()):
        cc=components(ss);hist[len(cc)]+=1;sizes[len(ss)]+=1
        if len(cc)>1:bad.append(dict(charge=list(key),states=sorted(ss),components=cc))
    return dict(margin=k,sources=[r['source'] for r in rows],states=sum(map(len,g.values())),nonempty_fibres=len(g),fibre_sizes=dict(sorted(sizes.items())),component_counts=dict(sorted(hist.items())),disconnected=bad)
def main():
    k=int(sys.argv[1]);r=census(k);p=HERE/f'census_k{k}.json'
    if '--write' in sys.argv:p.write_text(json.dumps(r,indent=2))
    elif json.loads(p.read_text())!=json.loads(json.dumps(r)):raise RuntimeError('census mismatch')
    print(json.dumps({k:v for k,v in r.items() if k not in ('sources','fibre_sizes','disconnected')}))
    if r['disconnected']:
        first=r['disconnected'][0];print('first hostile:',first['charge'],'component sizes',list(map(len,first['components'])))
if __name__=='__main__':main()
