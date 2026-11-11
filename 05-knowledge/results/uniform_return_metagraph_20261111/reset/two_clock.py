"""Exact two-module exponential response; probes splitting the growing prefix."""
from collections import defaultdict
from pathlib import Path
import json,math,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'chanpak_five_core_bridges_20261110/alternatives'))
import check as a
from classify_modules import H,strong,need

def multiruns(rows, modules=(0,1)):
    others=[v for v in range(len(rows)) if v not in modules]
    vertices=others+list(modules)
    R=a.induced(rows,vertices)
    n=len(others);full=(1<<n)-1
    level={(0,n+2,0,0):1};ans=defaultdict(int)
    # Between two ordinary vertices, a transitive pair supports at most2
    # different placeholders. There is no repeated same-module adjacency.
    for length in range(3*n+3):
        nxt=defaultdict(int)
        for (mask,last,r,s),count in level.items():
            if mask==full and r and s:ans[r,s]+=count
            for v in range(n+2):
                if v<n and mask>>v&1:continue
                if last!=n+2 and not(R[last]>>v&1):continue
                key=(mask|(1<<v) if v<n else mask,v,r+(v==n),s+(v==n+1))
                nxt[key]+=count
        level=nxt
    need(not level,'bounded multirun length')
    return ans

def exponential(runs,shifts=(0,0)):
    out=defaultdict(int)
    for (r,s),count in runs.items():
        for i in range(1,r+1):
            for j in range(1,s+1):
                out[i*j]+=count*(-1)**(r+s-i-j)*math.comb(r,i)*math.comb(s,j)*i**shifts[0]*j**shifts[1]
    return {k:v for k,v in sorted(out.items()) if v}

def inflate_pair(rows,m,n):
    labels=[(v,j) for v in range(len(rows)) for j in range(m if v==0 else n if v==1 else 1)]
    return tuple(sum((i!=j and ((v==w and p<q) or (v!=w and rows[v]>>w&1)))<<j
                     for j,(w,q) in enumerate(labels)) for i,(v,p) in enumerate(labels))

def main():
    old=HERE.parent.parent/'chanpak_five_core_bridges_20261110/alternatives/reset.json'
    records=json.loads(old.read_text())['replacement_probes'][1]['hits']
    probes=[]
    for rec in records:
        rows=rec['rows'];runs=multiruns(rows)
        e=exponential(runs)
        need(sum(k*c for k,c in e.items())==H(rows),'k1 two-module specialization')
        need(sum(c*k*k for k,c in e.items())==H(inflate_pair(rows,2,2)),'k2 literal two-module inflation')
        odd=exponential(runs,(0,1))
        need(sum(k*c for k,c in odd.items())==H(inflate_pair(rows,1,2)),'odd shift literal count')
        probes.append({'rows':rows,'strong':strong(rows),'outmasks':rec['outmasks'],
                       'even_exponential':e,'odd_exponential':odd,
                       'even_return':e=={1:5,4:34},'odd_return':odd=={1:5,4:68}})
    report={'universe':'All22 port-b-deletion/two-birth141 hits at original k2',
            'probes':probes,'strong_even_returns':sum(p['strong'] and p['even_return'] for p in probes),
            'strong_odd_returns':sum(p['strong'] and p['odd_return'] for p in probes)}
    (HERE/'two_clock.json').write_text(json.dumps(report,indent=2)+'\n')
    print('two-clock split',len(probes),'strong even',report['strong_even_returns'],'strong odd',report['strong_odd_returns'])
if __name__=='__main__':main()
