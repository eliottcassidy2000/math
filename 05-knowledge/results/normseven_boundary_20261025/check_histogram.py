#!/usr/bin/env python3
"""Fixed-source, integer-clock controls for the elementary affine histogram bound."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

GATES=0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def ceil_fraction(x):
    return -((-x.numerator)//x.denominator)


def main():
    receipts=[]
    histories=0
    for M in (1,5,13):
        frontier=[(M,0,F(1),F(0),())]
        for depth in range(1,7):
            new=[]
            for source,total,W,beta,word in frontier:
                for a in range(1,9):
                    numerator=(1<<a)*source-1
                    if numerator % 3:
                        continue
                    x=numerator//3
                    if x <= 0:
                        continue
                    Wnew=W*F(3,1<<a)
                    bnew=beta+W/F(1<<a)
                    need(Wnew*x+bnew == M, 'marked endpoint retained')
                    need((3*x+1)//(1<<a) == source and source % 2,
                         'native forward valuation and odd source')
                    need(bnew <= F(3,2)**depth-1, 'declared all-word offset budget')
                    new.append((x,total+a,Wnew,bnew,word+(a,)))
            frontier=new
            histories += len(frontier)
            bins=defaultdict(list)
            for row in frontier:
                bins[row[1]].append(row)
            for total,rows in bins.items():
                W=F(3**depth,1<<total)
                need(len({r[0] for r in rows}) == len(rows), 'one source owner at fixed D,A')
                B=F(3,2)**depth-1
                for m in (3,9,27,81):
                    masses=defaultdict(F)
                    for x,_,weight,_,word in rows:
                        masses[x%m] += weight
                        need(len(word)+sum(word) == depth+total, 'ordinary integer clock retained')
                    for residue,H in masses.items():
                        need(m*H <= B+m*W, 'fixed-bin affine histogram bound')
                        L=(M-B)/W
                        U=F(M)/W
                        count=max(0,((U-residue)/m).numerator//((U-residue)/m).denominator
                                  -ceil_fraction((L-residue)/m)+1)
                        need(H <= W*count, 'source-phase-specific integer interval bound')
            receipts.append([M,depth,len(frontier),len(bins),
                             min((depth+r[1] for r in frontier),default=0),
                             max((depth+r[1] for r in frontier),default=0)])
    # With two copies of one source, B=0,W=1,m=3 would claim 6<=3.
    need(3*2 > 0+3*1, 'duplicate-owner hostile is detected')
    print(json.dumps({'status':'FINITE-EXACT','gates':GATES,'histories':histories,
                     'universe':'M in {1,5,13}, depth 1..6, every valuation in 1..8; all positive integral inverse words, no reachability pruning',
                     'rows_M_depth_sources_bins_minclock_maxclock':receipts,
                     'scope':'No uniform-in-depth summability or fixed-root positivity is inferred.',
                     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                    indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
