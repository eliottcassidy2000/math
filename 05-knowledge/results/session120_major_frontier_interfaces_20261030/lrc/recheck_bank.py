#!/usr/bin/env python3
"""Recheck inherited finite-core bank, including its proved tail cutoffs."""
from pathlib import Path
import json
from math import ceil
from fractions import Fraction as F
from fibre_transfer import safe, transfer_accept, require

def run():
    bank=json.loads(Path(__file__).with_name('source_bank.json').read_text())
    r=F(1,14); n=0; result=[]
    for row in bank['cores']:
        c=tuple(row['core']); intervals=safe(c,r)
        W=max(b-a for a,b in intervals)
        require(W>0,'bank width')
        B=ceil(2*r/W)
        require(B==row['all_speeds_from'],'inherited tail cutoff')
        for v in range(1,B):
            require(transfer_accept(intervals,1,v,r),'base finite verdict')
            require(bool(safe(c+(v,),r)),'independent comb verdict')
            n+=1
        result.append({'core':list(c),'tail_cutoff':B,'finite_queries':B-1})
    require(n==1458 and len(result)==28,'declared bank universe')
    return {'status':'PASS','cores':result,'finite_queries':n,'infinite_tail':'For v >= ceil(2r/W), a width-W interval crosses the safe phase arc.', 'scope':'28 inherited core shapes; all dilations follow analytically, no new primitive shapes.'}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
