#!/usr/bin/env python3
"""Compare independent complete named cube with producer data and literal coarse controls."""
from pathlib import Path
from itertools import combinations
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
P=ROOT/'odd_pair_quotient';A=ROOT/'odd_pair_quotient_audit'

def need(ok,why):
    if not ok:raise RuntimeError(why)

def hp(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for S in range(1,1<<n):
        if S.bit_count()==1:continue
        for last in range(n):
            if not S>>last&1:continue
            before=S^(1<<last)
            d[S][last]=sum(d[before][v] for v in range(n) if before>>v&1 and a[v]>>last&1)
    return sum(d[-1])

def extend(a,x,c):
    n=len(a);rows=list(a)+[1<<(n+1),0]
    for i in range(n):
        if x>>i&1:rows[i]|=1<<n
        else:rows[n]|=1<<i
        if c>>i&1:rows[i]|=1<<(n+1)
        else:rows[n+1]|=1<<i
    return rows

def main():
    cube=json.loads((A/'full_named_cube.json').read_text());full=json.loads((P/'full_named6.json').read_text());basis=[a|(c<<6) for a,c in cube['basis']]
    source=next(s for s in cube['states'] if s['quotient_conversed']==0 and s['pair_bits']==0);target=next(s for s in cube['states'] if s['quotient_conversed']==1 and s['pair_bits']==0)
    need(source['rows']==full['T'] and target['rows']==full['U'],'physical rows')
    for state,key in [(source,'source_values'),(target,'target_values')]:
        v=state['vector'];need(len(v)==484,'complete physical coefficients')
        for m,expected in enumerate(full[key]):need(sum(x for s,x in zip(basis,v) if s&m==s)==expected,'all retained physical truth values')
    diff={s:x-y for s,x,y in zip(basis,source['vector'],target['vector']) if x!=y}
    need(diff=={r['mask']:r['coefficient'] for r in full['contrast_coefficients']},'all full named cubic/quartic coefficients')
    need(len({tuple(s['vector']) for s in cube['states']})==16,'exact 16-state cube')
    literal=0;records=[]
    for case in json.loads((P/'coarse_cases.json').read_text()):
        q=len(case['Q']);N=1<<q
        choices=sorted({0,N*N-1,min(N+1,N*N-1),min(3*N+2,N*N-1)})
        for m in choices:
            a,c=m&(N-1),m>>q;x=sum(3<<(2*i) for i in range(q) if a>>i&1);y=sum(3<<(2*i) for i in range(q) if c>>i&1)
            ft,fu=hp(extend(case['T'],x,y)),hp(extend(case['U'],x,y));need((ft,fu)==(case['actual_source_values'][m],case['actual_target_values'][m]),'independent coarse literal DP');literal+=2;records.append([q,m,ft,fu])
    producer=json.loads((P/'normal.stdout').read_text());need((P/'normal.stdout').read_bytes()==(P/'optimized.stdout').read_bytes(),'producer normal and optimized receipts')
    for file,digest in producer['source_sha256'].items():need(hashlib.sha256((P/file).read_bytes()).hexdigest()==digest,'producer source pin')
    for file,digest in producer['data_sha256'].items():need(hashlib.sha256((P/file).read_bytes()).hexdigest()==digest,'producer data pin')
    need(producer==json.loads((P/'receipt.json').read_text()),'saved producer receipt')
    return dict(status='PASS',physical_truth_values_compared=8192,complete_physical_coefficients_compared=968,complete_physical_contrast_nonzeros=len(diff),full_named_cube_states=16,independent_backward_dp_values=literal,literal_controls=records,source_sha256=producer['source_sha256'],data_sha256=producer['data_sha256'],primary_cycle_source='https://arxiv.org/pdf/2307.05569',primary_cycle_result='Theorem 1.39, arXiv v1, PDF page14; verified directly.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
