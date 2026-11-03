"""Minimize explicitly named partial re-markings; retain illegal commands."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
from scout import rows,h_count
HERE=Path(__file__).resolve().parent

def chart(n):
    edges=list(combinations(range(n),2));free=[(a,b) for a,b in edges if b>a+1]
    ef={e:k for k,e in enumerate(edges)};ff={e:k for k,e in enumerate(free)}
    fixed=sum(1<<ef[a,b] for a,b in edges if b==a+1)
    def expand(x):return fixed+sum(((x>>k)&1)<<ef[e] for k,e in enumerate(free))
    return free,ff,expand

def permutation_map(n,p,free,ff):
    # Compile guard and affine output from actual arc orientation.
    guard=val=0;impossible=False
    for a,b in zip(p,p[1:]):
        lo,hi=sorted((a,b))
        if hi==lo+1:
            if a>b:impossible=True
        else:
            guard|=1<<ff[lo,hi]
            if a<b:val|=1<<ff[lo,hi]
    output=[]
    for a,b in free:
        u,v=p[a],p[b];lo,hi=sorted((u,v))
        output.append((None,int(u<v)) if hi==lo+1 else (ff[lo,hi],int(u>v)))
    result=[]
    for x in range(1<<len(free)):
        if impossible or x&guard!=val:result.append(-1);continue
        result.append(sum((z if j is None else ((x>>j)&1)^z)<<k for k,(j,z) in enumerate(output)))
    return result

def minimize(outputs,maps):
    def compress(keys):
        names={};return [names.setdefault(k,len(names)) for k in keys]
    labels=compress([(h,tuple(t[x]>=0 for t in maps)) for x,h in enumerate(outputs)])
    rounds=[len(set(labels))]
    while True:
        nxt=compress([(labels[x],tuple(labels[t[x]] if t[x]>=0 else -1 for t in maps)) for x in range(len(outputs))])
        rounds.append(len(set(nxt)))
        if rounds[-1]==rounds[-2]:return nxt,rounds
        labels=nxt

def main():
    results=[]
    for n in range(3,8):
        free,ff,expand=chart(n);size=1<<len(free);top=size-1
        counts=[h_count(rows(n,expand(x))) for x in range(size)]
        C=[x^top for x in range(size)]
        rho=[sum(((x>>ff[n-1-b,n-1-a])&1)<<k for k,(a,b) in enumerate(free)) for x in range(size)]
        generators=[];names=[]
        for i in range(n):
            for j in range(i+2,n):
                p=list(range(n));v=p.pop(i);p.insert(j,v)
                t=permutation_map(n,p,free,ff)
                if any(y>=0 for y in t):generators.append(t);names.append(['right',i,j])
                inv=[-1]*size
                for x,y in enumerate(t):
                    if y>=0:inv[y]=x
                if any(y>=0 for y in inv):generators.append(inv);names.append(['left',j,i])
        row={'n':n,'cube':size,'relocation_labels':names}
        for label,maps in [('C_only',[C]),('C_relocations',[C]+generators),('C_relocations_converse',[C,rho]+generators)]:
            cls,rounds=minimize(counts,maps)
            row[label]={'rounds':rounds,'classes':len(set(cls)),'class_sizes':dict(sorted(Counter(Counter(cls).values()).items()))}
        cls,rounds=minimize([0]*size,[C]+generators)
        row['legality_only']={'rounds':rounds,'classes':len(set(cls))}
        if n==7:
            p=[0,5,6,3,4,1,2];R=permutation_map(n,p,free,ff)
            for label,maps in [('C_one_mark',[C,R]),('C_one_mark_converse',[C,R,rho])]:
                cls,rounds=minimize(counts,maps)
                row[label]={'rounds':rounds,'classes':len(set(cls)),
                            'class_sizes':dict(sorted(Counter(Counter(cls).values()).items()))}
            row['mark_witness']={'permutation':p,'source':24511,'target':R[24511]}
        results.append(row)
        print(json.dumps(row),flush=True)
    (HERE/'carrier_scout.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
