"""Exact strong unit-ready minimum: scalar, grammar, and literal controls."""
from pathlib import Path
import sys,json,hashlib
from math import comb
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'changing_core_returns_20261113'/'returns'))
import inverse as old

def minimum(M):
    r=M-4
    if r<0:raise ValueError('strong unit-ready order must be at least4')
    if r==0:return 5,[1]
    if r==1:return 9,[2]
    if r==2:return 17,[4]
    if r==3:return 33,[8,9]
    s=(r-1)//2;n=(1<<r)+(1<<s)-1
    return 3*(1<<r)+3*(1<<s)+(1<<(r-s))-1,[n]

def scalar(n):return 3*n+old.g(n)+(n&1)

def main():
    scalar_rows=[]
    for r in range(19):
        values=[scalar(n) for n in range(1<<r,1<<(r+1))]
        best=min(values);addresses=[(1<<r)+i for i,x in enumerate(values) if x==best]
        old.need((best,addresses)==minimum(r+4),'complete scalar minimum and equality addresses')
        scalar_rows.append(dict(order=r+4,H=best,addresses=addresses))
    census=RESULTS.parent.parent/'session212/cycle-spectrum/tournaments8.txt'
    strong_controls=0;triangle_min=None;p2_min=None
    for line in census.read_text().splitlines():
        T=tuple(map(int,line.split()))
        if not old.f.a.strong(T):continue
        n=len(T);scores=sorted(row.bit_count() for row in T)
        old.need(all(sum(scores[:k])>=comb(k,2)+1 for k in range(1,n)),'strong score-prefix gate')
        triangles=comb(n,3)-sum(comb(x,2) for x in scores)
        old.need(triangles>=n-2,'triangle minimum')
        H,_=old.f.a.tables(T);full=(1<<n)-1
        P2=sum(H[mask]*H[full^mask] for mask in range(full+1))
        old.need(P2>=(triangles+2)*(1<<(n-1))>=n*(1<<(n-1)),'OCF two-bin floor')
        strong_controls+=1;triangle_min=triangles if triangle_min is None else min(triangle_min,triangles);p2_min=P2 if p2_min is None else min(p2_min,P2)
        old.f.a.tables.cache_clear()
    grammar=[];total=0
    for n in range(3,16):
        best=None;min_ranks=[]
        for rank in range(old.f.a.class_counts(n)[n]):
            c=old.f.compile_rank(n,rank);value=c[0]+c[4]-1;total+=1
            if best is None or value<best:best=value;min_ranks=[rank]
            elif value==best:min_ranks.append(rank)
        expected,addresses=minimum(n+1)
        old.need(best==expected and len(min_ranks)==len(addresses),'complete independent grammar minimum')
        grammar.append(dict(order=n+1,H=best,parent_ranks=min_ranks))
    literals=[];single=old.cores(1)[0]
    for M in range(4,16):
        H,addresses=minimum(M)
        for n in addresses:
            T,y,V,x,w=old.physical(single,n)
            old.need(len(T)==M and old.f.a.strong(T),'literal extremal size/strongness')
            old.need(old.direct(T)==H,'literal extremal full H')
            old.need(old.direct(V,x)==1 and old.direct(V)-old.direct(T,y)==1,'literal extremal readiness')
            literals.append(dict(order=M,H=H,address=n,word=''.join(map(str,w)),rows=T,root=y,predecessor=x,degree_multiset=sorted(row.bit_count() for row in T)))
    rec=dict(status='PASS',scalar_exhaustive_max_r=18,scalar_minima=scalar_rows,strong_order8_controls=strong_controls,minimum_triangle_count=triangle_min,minimum_P2=p2_min,source_census_sha256=hashlib.sha256(census.read_bytes()).hexdigest(),complete_grammar_parent_order=15,complete_grammar_classes=total,grammar_minima=grammar,literal_extremizers=literals)
    suffix='' if __debug__ else '.optimized'
    (HERE/('unit_ready_minimum_receipt'+suffix+'.json')).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS exact strong unit-ready minimum:',strong_controls,'strong literal contexts;',total,'grammar classes')

if __name__=='__main__':main()
