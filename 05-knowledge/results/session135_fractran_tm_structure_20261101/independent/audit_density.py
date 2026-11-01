#!/usr/bin/env python3
"""Exact cylinder weights using independently detected PRIMEGAME basins."""
import itertools,json,math,sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
from audit_early import independent_basin,check,P
ROOT=Path(__file__).resolve().parents[1]/'primegame_machine'
def exact(obj):return F(obj['numerator'],obj['denominator'])
def main():
    cap=1;d=math.prod(p**(cap+1) for p in P);counts=dict.fromkeys(('prime','cycle2','cycle6'),0);nums=dict.fromkeys(counts,0)
    for x in itertools.product(range(cap+1),repeat=10):
        label=independent_basin(x);counts[label]+=1;nums[label]+=math.prod((p-1)*p**(cap-e) for p,e in zip(P,x))
    saved=json.loads((ROOT/'APPROXIMATION.normal.json').read_text());tail=d-math.prod(p**(cap+1)-1 for p in P)
    check(saved['common_denominator']==d and saved['error_numerator']==tail and saved['census']==counts,'density cap certificate')
    for key in counts:
        check(saved['intervals'][key]=={'lower_numerator':nums[key],'upper_numerator':nums[key]+tail},('density sum',key))
    controls=(11,13,17,19,23,29);c0=math.prod(F(p-1,p) for p in controls)
    onehot=c0*(1+sum(F(1,p) for p in controls))
    prime=c0*(F(2,3)+(F(1,11)+F(1,13))*F(9,35)+(F(1,17)+F(1,19)+F(1,29))*F(3,5)+F(1,23))
    extra=F(2,5)*(F(22,23)-c0*(1+sum(F(1,p) for p in controls if p!=23)))
    six=F(18944,46835)*sum(F(1,n) for n in (4301,17765,16302,9867,8671,23023))
    lower={'prime':prime,'cycle2':onehot-prime+extra,'cycle6':six};report=json.loads((ROOT/'DENSITY.normal.json').read_text())
    for key,value in lower.items():
        check(exact(report['basin_intervals'][key]['lower'])==value,('analytic lower',key))
        check(exact(report['basin_intervals'][key]['upper'])==1-sum(v for k,v in lower.items() if k!=key),('analytic upper',key))
    for m in range(65):
        H=m+4;error=1-math.prod(1-F(1,p**(H+1)) for p in P)
        check(error<F(1,2**m),'precision tail')
    print(json.dumps({'status':'PASS','independent_cap':cap,'states':2**10,'census':counts,'common_denominator':d,'analytic_lower_bounds':{k:str(v) for k,v in lower.items()},'precision_tail_controls':'m=0..64'},indent=2))
if __name__=='__main__':main()
