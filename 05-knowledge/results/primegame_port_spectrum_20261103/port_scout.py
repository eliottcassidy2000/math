"""Exact one-port birth and two-port commutator probes; standard library."""
from itertools import combinations,permutations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def rows(n,mask):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if mask>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return a

def encode(a):return sum(((a[i]>>j)&1)<<k for k,(i,j) in enumerate(combinations(range(len(a)),2)))

def table(a):
    n=len(a);end=1<<n;starts=[[0]*n for _ in range(end)];hs=[0]*end;hs[0]=1
    for s in range(1,end):
        for v in range(n):
            if not s>>v&1:continue
            rest=s^(1<<v)
            if not rest:starts[s][v]=1
            else:starts[s][v]=sum(starts[rest][w] for w in range(n) if rest>>w&1 and a[v]>>w&1)
        hs[s]=sum(starts[s])
    return hs,starts

def h(a):return table(a)[0][-1]
def add(a,x):
    n=len(a)
    return [r|((1<<n) if i!=x else 0) for i,r in enumerate(a)]+[1<<x]
def converse(a):return [((1<<len(a))-1)^(1<<i)^r for i,r in enumerate(a)]

def classes(n):
    unseen=set(range(1<<(n*(n-1)//2)));ps=list(permutations(range(n)));result=[]
    while unseen:
        m=min(unseen);a=rows(n,m);orbit=set()
        for p in ps:orbit.add(encode([sum(((a[p[i]]>>p[j])&1)<<j for j in range(n)) for i in range(n)]))
        unseen.difference_update(orbit);result.append((m,a,len(orbit)))
    return result

def port_counts(a):
    n=len(a);full=(1<<n)-1;hs,starts=table(a);_,ends=table(converse(a))
    return [hs[full]-ends[full][x]+sum(starts[s][x]*hs[full^s] for s in range(1,full+1) if s>>x&1) for x in range(n)]

def main():
    report=[]
    for n in range(2,7):
        entries=[]
        for m,a,orbit in classes(n):
            H=h(a);bs=port_counts(a)
            if bs!=[h(add(a,x)) for x in range(n)]:raise RuntimeError('birth formula')
            entries.append({'mask':m,'orbit':orbit,'H':H,'children':bs,'binary_ports':{str(d):[i for i,b in enumerate(bs) if b==2*H+d] for d in (-1,1)}})
        report.append({'n':n,'classes':len(entries),'entries':entries})
    (HERE/'port_scout.json').write_text(json.dumps(report,indent=2)+'\n')
    for r in report:
        es=r['entries'];both=[e for e in es if all(e['binary_ports'].values())]
        print(r['n'],len(es),'both',len(both),'H',sorted({e['H'] for e in both}))

if __name__=='__main__':main()
