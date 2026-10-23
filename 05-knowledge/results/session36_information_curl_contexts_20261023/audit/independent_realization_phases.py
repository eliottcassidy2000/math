#!/usr/bin/env python3
"""Independent Boolean bitsets and GF2 elimination, no parity-DSU import."""
from independent_audit import mat,need
import itertools,json
from collections import defaultdict,Counter
from pathlib import Path
OUT=Path(__file__).resolve().parent

def c3mask(z):
    n=len(z);return sum(all(sum(z[v][w]==1 for w in s)==1 for v in s)<<b for b,s in enumerate(itertools.combinations(range(n),3)))

def positive_rank(n,target):
    pairs=list(itertools.combinations(range(n),2));ix={p:i for i,p in enumerate(pairs)};m=len(pairs);rows=[]
    for b,(i,j,k) in enumerate(itertools.combinations(range(n),3)):
        if target>>b&1:
            rows.extend([(1<<ix[i,j])^(1<<ix[j,k]),(1<<ix[i,j])^(1<<ix[i,k])^(1<<m)])
    rank=0
    for col in range(m):
        pivot=next((i for i in range(rank,len(rows)) if rows[i]>>col&1),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        for i in range(len(rows)):
            if i!=rank and rows[i]>>col&1:rows[i]^=rows[rank]
        rank+=1
    inconsistent=any(r==(1<<m) for r in rows)
    return inconsistent,m-rank

def main():
    rows=[];specific={}
    for n in range(1,6):
        m=n*(n-1)//2;triples=list(itertools.combinations(range(n),3));groups=defaultdict(list);positive=[0]*len(triples)
        for code in range(1<<m):
            h=c3mask(mat(n,code));groups[h].append(code)
            for b in range(len(triples)):
                if h>>b&1:positive[b]|=1<<code
        counts=Counter();phasechecks=0;phaseconsistent=Counter();phasereal=Counter();phaseprime=Counter()
        for h in range(1<<len(triples)):
            feasible=(1<<(1<<m))-1
            for b in range(len(triples)):
                if h>>b&1:feasible &= positive[b]
            bad,free=positive_rank(n,h)
            need(bool(feasible)==(not bad),('GF2/positive enumeration',n,h))
            if feasible:need(feasible.bit_count()==1<<free,('phase dimensions',n,h,free))
            state='XOR_inconsistent' if bad else ('realizable' if h in groups else 'residual_NAE_unsatisfiable')
            counts[state]+=1
            if not bad:phaseconsistent[free]+=1
            if h in groups:
                phasereal[free]+=1
                need(all((code^((1<<m)-1)) in groups[h] for code in groups[h]),('converse pairing',n,h))
                z=mat(n,groups[h][0]);prime=True
                for sm in range(1,1<<n):
                    inside=[j for j in range(n) if sm>>j&1];outside=[j for j in range(n) if not sm>>j&1]
                    if 1<len(inside)<n and all(len({z[v][w] for v in inside})==1 for w in outside):prime=False;break
                if prime:
                    phaseprime[free]+=1
                    need(len(groups[h])==(1 if n==1 else 2),('indecomposable fibre',n,h))
            if n==5 and h in(13,12):specific[h]=dict(positive_triples=[s for b,s in enumerate(triples) if h>>b&1],free_phases=free,solutions=groups.get(h,[]))
            # Literal model: each directed cyclic triple has literals (y_ij,y_jk,not y_ik).
            for code in groups.get(h,[]):
                ix={p:i for i,p in enumerate(itertools.combinations(range(n),2))}
                for b,(i,j,k) in enumerate(triples):
                    lits=[code>>ix[i,j]&1,code>>ix[j,k]&1,1-(code>>ix[i,k]&1)]
                    need((len(set(lits))==1)==bool(h>>b&1),('literal convention',n,h,code));phasechecks+=1
        row=dict(n=n,hypergraphs=1<<len(triples),counts=dict(counts),literal_checks=phasechecks,phase_histogram_consistent=dict(phaseconsistent),phase_histogram_realizable=dict(phasereal),indecomposable_phase_histogram=dict(phaseprime));rows.append(row);print(row,flush=True)
    need(specific[12]['positive_triples']==[(0,1,4),(0,2,3)],'wrong minimal hostile mask')
    need(specific[12]['free_phases']==6 and specific[12]['solutions']==[],'wrong hostile phases')
    need(specific[13]['free_phases']==4 and specific[13]['solutions']==[117,906],'wrong indecomposable phases')
    report=dict(status='PASS',rows=rows,witnesses=specific,method='Literal tournament bitsets independently cross-checked by GF2 Gaussian elimination; no DSU or producer imports.')
    (OUT/'realization_phase_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
