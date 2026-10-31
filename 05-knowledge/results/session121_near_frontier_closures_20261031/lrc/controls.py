"""Closed endpoints, hostile joint incompatibility, and conductor equality controls."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import ceil
from pathlib import Path
import json
import source_lrc_banks as literal
from triple_certificate import attach

def run():
    checks=0
    def need(ok,tag):
        nonlocal checks
        checks+=1
        if not ok: raise RuntimeError(tag)
    for r in (F(1,14),F(1,6),F(1,4),F(1,3),F(1,2)):
        for a in range(1,9):
            I=literal.safe_set((a,),r)
            for b in range(1,13):
                need(attach(I,b,r)==literal.safe_set((a,b),r),'complete clipped interval vs full comb')
    I=literal.safe_set((1,2),F(1,4))
    need(bool(attach(I,3,F(1,4))) and bool(attach(I,12,F(1,4))),'hostile individual attachments')
    need(not attach(attach(I,3,F(1,4)),12,F(1,4)),'hostile joint rejection')
    need(attach([(F(1,3),F(1,3))],1,F(1,3))==[(F(1,3),F(1,3))],'closed singleton preserved')
    need(attach([(F(1,3),F(1,3))],3,F(1,3))==[],'unsafe singleton rejected')
    conductor=0
    for k in (2,3,4):
        r=F(1,4*k-2)
        for W in (F(1,2),F(1,3),F(2,5)):
            B=ceil(1/W)
            for left in (F(0),F(1,17),F(2,7)):
                for speeds in combinations_with_replacement(range(B,B+6),k):
                    I=[(left,left+W)]
                    for v in speeds: I=attach(I,v,r)
                    need(sum(h-l for l,h in I)>0,'conductor at exact radius threshold')
                    conductor+=1
    return {'checks':checks,'complete_comb_checks':480,'conductor_boundary_controls':conductor,'hostile_and_point_controls':4}
if __name__=='__main__':
    print(json.dumps(run(),sort_keys=True))
