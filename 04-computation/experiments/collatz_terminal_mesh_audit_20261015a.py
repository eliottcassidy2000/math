#!/usr/bin/env python3
"""Independent literal-source and formal-terminal checks for the mean ledger."""
import itertools
import json
from fractions import Fraction as F
from pathlib import Path
from collatz_terminal_mesh_20261015a import certificate,endpoints,T
from collatz_clock_mesh_audit_20261015a import need

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'05-knowledge/results'

def formal_checks():
    count=0
    for depth in (0,1,2,3,8):
        k=2;X=2**k;N=2*X
        row=certificate(k,depth,keep_ledger=True);ends,_=endpoints(N,depth)
        keys=sorted({e for e,s in ends[1:] if e!=1})
        for vals in itertools.product((-1,0,1),repeat=len(keys)):
            z=dict(zip(keys,vals));z[1]=1
            t=[0]+[s*z[e] for e,s in ends[1:]]
            J=0;currents=[]
            for n in range(1,N):
                if n>=X:currents.append(J)
                if n+1<N:J+=F((0,-1,1)[n%3]*(t[n+1]-t[n]),2)
            b=row['b']
            for i,coeff in enumerate(row['terminal_weight_ledger']):
                actual=2*sum(currents[i*b:(i+1)*b])
                need(actual==sum(v*z[int(e)] for e,v in coeff.items()),'formal terminal row versus direct current')
                known,mass=row['block_known_twice_mean_and_residual_l1'][i]
                need(abs(actual-known)<=mass,'formal terminal bound')
                count+=1
    return count

def main():
    rows=[]
    for name in ('collatz_terminal_mesh_20261015a.json','collatz_terminal_mesh_large_20261015a.json'):
        rows+=json.loads((OUT/name).read_text())['rows']
    K=max(r['k'] for r in rows);N=2**(K+1)
    # Different stopping condition: first descent to an earlier source,
    # never the fixed-depth terminal controller used to construct bounds.
    p=bytearray(N+1)
    for n in range(2,N+1):
        x=T(n);steps=1
        while x>=n:
            x=T(x);steps+=1;need(steps<1000000,'independent descent cap')
        p[n]=(steps+p[x])%2
    J=0;actual={};block_sum=0
    for n in range(1,N):
        k=n.bit_length()-1;X=2**k;b=2**(k//2)
        block_sum+=J
        if (n-X+1)%b==0:
            actual.setdefault(k,[]).append(block_sum);block_sum=0
        J+=(0,-1,1)[n%3]*((1-2*p[n+1])-(1-2*p[n]))//2
    native={r['k']:r for r in json.loads((OUT/'collatz_clock_mesh_20261015a.json').read_text())['octaves']}
    prev={};checks=0
    for row in rows:
        k=row['k'];b=row['b'];means=actual[k]
        need(len(means)==len(row['block_known_twice_mean_and_residual_l1']),'physical mesh size')
        coarse=sum(F(m*m,b) for m in means)
        nr=native[k]['meshes'][0]
        need(coarse==F(int(nr['coarse_numerator']),nr['coarse_denominator']),'native coarse control')
        B=[]
        for m,(known,mass) in zip(means,row['block_known_twice_mean_and_residual_l1']):
            need(abs(2*m-known)<=mass,'actual marked source within terminal enclosure')
            B.append(abs(known)+mass);checks+=1
        need(sum(F(v*v,4*b) for v in B)==F(*row['coarse_upper_exact']),'upper sum exact')
        if k in prev:
            need(all(v<=old for v,old in zip(B,prev[k])),'controller monotonicity at every physical block')
        prev[k]=B
    result={'status':'FINITE-EXACT independent audit','literal_descent_sources':N,
            'actual_block_enclosures_checked':checks,'formal_terminal_assignment_checks':formal_checks(),
            'controller_rows':len(rows)}
    (OUT/'collatz_terminal_mesh_audit_20261015a.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
