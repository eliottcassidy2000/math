#!/usr/bin/env python3
"""Source-preserving signed-terminal upper bounds for coarse current means.
No unknown terminal is assigned a guessed sign. Everything is integer/rational.
"""
import json
import sys
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from collatz_clock_mesh_audit_20261015a import envelope,need

def T(n):return (3*n+1)//2 if n%2 else n//2

def endpoints(N,depth):
    result=[None]*N;unresolved=0
    for n in range(1,N):
        x=n;s=1
        for j in range(depth):
            if x==1:break
            x=T(x);s=-s
        result[n]=(x,s)
        unresolved+=x!=1
    return result,unresolved

def certificate(k,depth,j=0,keep_ledger=False):
    X=2**k;b=2**max(0,k//2-j);N=2*X
    ends,unresolved=endpoints(N,depth)
    prefix=defaultdict(int);next_source=1;C=F(0);R=F(0)
    terminal_union=set();blocks=[];ledger=[]
    for a in range(X,2*X,b):
        while next_source<a:
            n=next_source;e,s=ends[n];prefix[e]+=(2 if n%3==2 else -1)*s;next_source+=1
        coeff={e:-b*v for e,v in prefix.items() if v}
        for n in range(a,a+b):
            e,s=ends[n]
            w=2 if n%3==2 else -1
            c=(-w*(a+b-n)+(0,-1,1)[n%3])*s
            coeff[e]=coeff.get(e,0)+c
        known=coeff.get(1,0);mass=sum(abs(v) for e,v in coeff.items() if e!=1)
        terminal_union.update(e for e,v in coeff.items() if e!=1 and v)
        C+=F((abs(known)+mass)**2,4*b);R+=envelope(a,b)
        blocks.append([known,mass])
        if keep_ledger:ledger.append({str(e):v for e,v in coeff.items() if v})
    Z=C+R;A=2*X*X-4*Z-X
    result={'k':k,'depth':depth,'j':j,'b':b,'source_count':N-1,
            'sources_not_grounded_by_this_controller':unresolved,
            'nonzero_unresolved_terminals':len(terminal_union),
            'coarse_upper_exact':[C.numerator,C.denominator],
            'coarse_upper_over_X2':float(C/(X*X)),
            'current_energy_upper_over_X2':float(Z/(X*X)),
            'certifies_E_le_2X2':A>=0 and A*A>=16*X*Z,
            'block_known_twice_mean_and_residual_l1':blocks}
    if keep_ledger:result['terminal_weight_ledger']=ledger
    return result

def main():
    rows=[]
    large=len(sys.argv)>2 and sys.argv[2]=='--large'
    for k in ((18,20) if large else (10,12,14,16)):
        for depth in ((160,192,224,256) if large else (16,32,64,96,128,160,192,256)):
            row=certificate(k,depth);rows.append(row)
            print({key:row[key] for key in ('k','depth','sources_not_grounded_by_this_controller','current_energy_upper_over_X2','certifies_E_le_2X2')},flush=True)
            if large and row['certifies_E_le_2X2']:break
    Path(sys.argv[1]).write_text(json.dumps({'status':'FINITE-EXACT source-weight certificates with unresolved terminals retained','rows':rows},indent=2)+'\n')

if __name__=='__main__':main()
