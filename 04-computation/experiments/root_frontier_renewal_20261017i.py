#!/usr/bin/env python3
"""Finite-horizon frontier renewal with literal source and shortcut-clock checks.

No loop waits for eventual entry. Every comparison is capped at 12 native
stages, and the moving remainder is retained. All arithmetic is integer.
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path


def need(p, msg):
    if not p: raise RuntimeError(msg)


def primitive(x):
    j=0
    while x%8==5:
        x=(x-1)//4;j+=1
    return x,j


def move(s):
    a=2 if s%8==1 else 1
    t,j=primitive((3*s+1)//2**a)
    return t,(-1)**a,a,j


def clean(v): return {s:c for s,c in v.items() if c}


def push(v, cutoff=None):
    out=defaultdict(int)
    for s,c in v.items():
        if cutoff is not None and s<=cutoff:out[s]+=c
        else:
            t,sg,_,_=move(s);out[t]+=sg*c
    return clean(out)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    horizon=12; renewal_checks=stopped_checks=clock_checks=energy_checks=0
    for Y in (1,3,10,100):
        low=[s for s in range(1,Y+1,2) if s%8!=5]
        returns=[{} for _ in range(horizon+1)]
        for s in low:
            t=s;sg=1
            for ell in range(1,horizon+1):
                t,a,_,_=move(t);sg*=a
                if t<=Y:
                    returns[ell][s]=(t,sg);break
        for N in range(1,65):
            sources=[];initial=defaultdict(int)
            for m in range(1,4*N+1,2):
                c=-3*(0,1,-1)[m%3] if m<2*N else (2 if m%3==2 else -1)
                s,j=primitive(m);initial[s]+=c
                # primitive seed, weight, literal endpoint, literal clock,
                # homogeneous clock, initial braid, terminal braid.
                sources.append([s,c,m,0,0,j,j])
            initial=clean(initial);full=initial;stopped=initial
            high={s:c for s,c in initial.items() if s>Y}
            frozen={s:c for s,c in initial.items() if s<=Y}
            low_history=[];E0=sum(c*c for c in initial.values());cross=0
            for d in range(horizon+1):
                actual_low={s:c for s,c in full.items() if s<=Y}
                arrival=frozen if d==0 else newly_arrived
                rhs=defaultdict(int,arrival)
                for ell in range(1,d+1):
                    for s,c in low_history[d-ell].items():
                        if s in returns[ell]:
                            t,sg=returns[ell][s];rhs[t]+=sg*c
                need(clean(rhs)==actual_low,'first-return renewal identity')
                renewal_checks+=1;low_history.append(actual_low)
                direct=defaultdict(int)
                for s,c,*_ in sources:direct[s]+=c
                need(clean(direct)==stopped,'source-specific stopped aggregation')
                stopped_checks+=1
                need(sum(c*c for c in stopped.values())==E0+2*cross,'finite-horizon energy ledger')
                energy_checks+=1
                if d==horizon:break
                full=push(full)
                moving={s:c for s,c in stopped.items() if s>Y}
                frozen={s:c for s,c in stopped.items() if s<=Y}
                land=defaultdict(list)
                for s,c in moving.items():
                    t,sg,_,_=move(s);land[t].append(sg*c)
                for t,cs in land.items():
                    cross+=(sum(cs)**2-sum(c*c for c in cs))//2
                    cross+=frozen.get(t,0)*sum(cs)
                stopped=push(stopped,Y)
                high_next=push(high)
                newly_arrived={s:c for s,c in high_next.items() if s<=Y}
                high={s:c for s,c in high_next.items() if s>Y}
                for row in sources:
                    s,c,x,tau,A,j0,jold=row
                    if s<=Y:continue
                    t,sg,a,j=move(s)
                    xnext=3*x+1;duration=0
                    while xnext%2==0:xnext//=2;duration+=1
                    need(duration==a+2*jold,'native-to-literal clock increment')
                    tau+=duration;A+=a+2*j
                    need(primitive(xnext)==(t,j),'literal endpoint braid')
                    need(tau==A+2*j0-2*j,'integer clock boundary correction')
                    row[:]=[t,c*sg,xnext,tau,A,j0,j];clock_checks+=1
    closed_frontiers=[]
    for Y in (3,10,100,1000):
        rank={1:0}
        for s in range(1,Y+1,2):
            if s%8==5:continue
            path=[];on_path=set();t=s
            while t not in rank:
                need(t not in on_path,'finite nucleus cycle control')
                need(len(path)<10000,'finite nucleus horizon cap')
                path.append(t);on_path.add(t);t=move(t)[0]
            while path:
                t=path.pop();rank[t]=1+rank[move(t)[0]]
        need(all(move(t)[0] in rank for t in rank),'forward-closed frontier')
        need(all(rank[move(t)[0]]<r for t,r in rank.items() if t!=1),'finite convergence ranks')
        closed_frontiers.append(dict(Y=Y,seeds=len(rank),max_seed=max(rank),max_rank=max(rank.values()),
                                     vertices=[[t,r] for t,r in sorted(rank.items())]))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    result=dict(status='FINITE-EXACT',sources='N=1..64',frontiers=[1,3,10,100],horizon=horizon,
                renewal_checks=renewal_checks,source_aggregation_checks=stopped_checks,
                finite_energy_checks=energy_checks,literal_clock_checks=clock_checks,
                closed_frontiers=closed_frontiers,
                prefix_hostile={'N':1,'Y_min':3,'coefficients':[[1,-3],[3,-1]],'energy':10},
                scope='finite horizons; moving complement retained; no universal entry assertion')
    args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
