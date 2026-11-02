"""Reproduce: python3 audit_tournaments.py > normal.json (also python3 -O).

Finite universe all 33,867 labelled nonempty tournaments through order six;
the n=6 unit graph has all strong children, all strong deletion roots,
single cut flips and opposite-bit two-coordinate exchanges, edges |dH|<=2.
Direct odd-cycle packing checks through n=5; source/sink, transitive,
three-cycle, and disconnected unit-graph controls are explicit.
"""
import json
from collections import Counter,defaultdict
from itertools import combinations
from tournament_kernel import *

def unit_components(n, hs, strong_codes):
    pairs=list(combinations(range(n),2));idx={p:k for k,p in enumerate(pairs)}
    parent={c:c for c in strong_codes}
    def find(c):
        while parent[c]!=c:
            parent[c]=parent[parent[c]];c=parent[c]
        return c
    def join(c,d):
        if d in parent and abs(hs[c]-hs[d])<=2:
            x,y=find(c),find(d)
            if x!=y:parent[max(x,y)]=min(x,y)
    edges=0
    for c in strong_codes:
        a=matrix(n,c)
        for r in range(n):
            if not strong(induced(a,[i for i in range(n) if i!=r])):continue
            for v in range(n):
                if v==r:continue
                d=c^(1<<idx[tuple(sorted((r,v)))])
                join(c,d)
            for v,w in combinations([i for i in range(n) if i!=r],2):
                if a[r][v]!=a[r][w]:
                    d=c^(1<<idx[tuple(sorted((r,v)))])^(1<<idx[tuple(sorted((r,w)))])
                    join(c,d)
    groups=defaultdict(list)
    for c in strong_codes:groups[find(c)].append(c)
    return sorted([{'size':len(cs),'H':sorted({hs[c] for c in cs}),
                    'representative_code':min(cs)} for cs in groups.values()],key=lambda d:d['H'])

def main():
    tables={};summaries=[];cycles_checked=0;literal_checked=0
    for n in range(1,7):
        hs=[];ss=[]
        for c in range(1<<(n*(n-1)//2)):
            a=matrix(n,c);h=hp(a);hs.append(h)
            require(h%2==1,'Redei parity')
            if strong(a):ss.append(c)
            if n<=5:
                require(h==literal_hp(a),'DP/literal mismatch');literal_checked+=1
                p=cycle_packing_polynomial(odd_cycles(a))
                require(h==sum(x*2**k for k,x in enumerate(p)),'odd-cycle formula');cycles_checked+=1
        tables[n]=(hs,ss)
        summaries.append({'n':n,'labelled':len(hs),'strong':len(ss),'H':sorted(set(hs))})
    require([x['strong'] for x in summaries]==[1,0,2,24,544,22320],'strong count')
    gaps=0;nested=0
    for n in range(1,5):
        for c in range(1<<(n*(n-1)//2)):
            b=matrix(n,c);s,e,q=exposed(b)
            const=[s[i]+e[i]+sum(q[j][i] for j in range(n)) for i in range(n)]
            cap=lambda cut:[const[i]+sum((q[i][j]-q[j][i])*((cut>>j)&1) for j in range(n)) for i in range(n)]
            for A in range(1<<n):
                a=append(b,A);aa=cap(A)
                for C in range(1<<n):
                    # x->y means new y loses to x, so y's x bit=0.
                    xy=append(a,C);yx=append(a,C|(1<<n))
                    hxy=tables[n+2][0][encoding(xy)]
                    hyx=tables[n+2][0][encoding(yx)]
                    pred=sum((((C>>i)&1)-((A>>i)&1))*aa[i] for i in range(n))
                    require(hxy-hyx==pred,'inter-ear orientation gap')
                    require(pred==sum((((C>>i)&1)-((A>>i)&1))*cap(C)[i] for i in range(n)),'skew affine identity')
                    gaps+=1
                    if A!=C and A&C==A:
                        floor=2*sum(hp(induced(b,[j for j in range(n) if j!=i])) for i in range(n) if (C^A)>>i&1)
                        require(pred>=floor,'nested positive-capacity floor');nested+=1
    comps=unit_components(6,*tables[6])
    require([x['H'] for x in comps]==[[15,17,19],[23,25],[27,29,31,33],[37],[41,43,45]],'unit components')
    controls={'transitive_H':hp(matrix(6,0)),'C3_H':hp(matrix(3,2)),
              'C3_to_C3_H':9,'unit_25_27_separated':True}
    c3=matrix(3,2);joined=tuple(row+(1,1,1) for row in c3)+tuple((0,0,0)+row for row in c3)
    require(hp(joined)==9,'SCC product control')
    print(json.dumps({'scope':'FINITE-EXACT; no all-order spectrum conclusion',
                      'labelled_summary':summaries,'literal_crosschecks':literal_checked,
                      'cycle_packing_crosschecks':cycles_checked,'interear_gap_cases':gaps,
                      'nested_positive_floor_cases':nested,'order6_unit_components':comps,
                      'controls':controls},sort_keys=True,indent=2))
if __name__=='__main__':main()
