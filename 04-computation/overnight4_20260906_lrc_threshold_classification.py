#!/usr/bin/env python3
"""Exact finite bases for the all-height old-threshold classification.

The full table was independently matched by literal intervals; this program
independently enumerates the two small coefficient-identity bases. The
all-height low/generic reductions are proved in the linked native report.
"""
import csv
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path
import sys

GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def main():
    sys.stdout.reconfigure(newline='\n')
    p=Path(__file__).resolve().parents[1]/'05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv'
    need(sha256(p.read_bytes()).hexdigest()=='c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa','pinned independently audited head')
    table={}
    for row in csv.DictReader(p.open(encoding='utf-8'),delimiter='\t'):
        w=tuple(int(row[x]) for x in ('a','b','c'))
        D=int(row['denominator'])
        E=tuple(Q(int(row['E'+str(i)+'_numerator']),D) for i in range(3))
        table[w]=(Q(int(row['mass_numerator']),D),min(E))
    need(len(table)==10074,'complete table cardinality')
    target=Q(6,77)
    additive=[]
    for c in range(1,62):
        for a in range(1,(c+1)//2):
            b=c-a
            if a>=b or a%3==0 or b%3==0 or c%3==0 or gcd(a,b)!=1:continue
            w=(a,b,c);additive.append(w)
            expected=Q(1,28) if w==(1,4,5) else None
            for value in table[w]:need(value==expected if expected is not None else value>target,'additive complete finite base '+str(w))
    need(len(additive)==146,'additive head61 size')
    norm4=[]
    units=[x for x in range(1,35) if x%3]
    for w in combinations(units,3):
        a,b,c=w
        if gcd(gcd(a,b),c)!=1 or not(c==2*a+b or c==a+2*b or 2*b==a+c):continue
        norm4.append(w)
        for value in table[w]:
            if w==(2,11,20):need(value==Q(11,140)>target,'norm4 unique above witness')
            elif w==(1,5,11):need(value==target,'norm4 unique equality witness')
            else:need(value<target,'norm4 strict finite remainder '+str(w))
    need(len(norm4)==88,'norm4 head34 size')
    need(Q(3,49)+Q(4,7*35)==Q(19,245)<target,'strict norm4 tail from35')
    need(Q(9,98)-Q(6,7*62)==Q(237,3038)>target,'strict additive tail from62')
    need(Q(46,665)<target,'incoming norm5 sharp bound is strict')
    print('FINITE-EXACT inherited literal head subbases; all-height statements require the analytic proof')
    print('additive_head61=146; both values above6/77 except(1,4,5), both1/28')
    print('norm4_head34=88; both above only(2,11,20), both equal only(1,5,11)')
    print('strict_tail_margins norm4_c35='+str(target-Q(19,245))+' additive_c62='+str(Q(237,3038)-target))
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)

if __name__=='__main__':main()
