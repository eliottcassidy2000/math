#!/usr/bin/env python3
"""Exact quantized fee transfer and all remaining E12 defect partitions."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
import json


def need(ok,message):
    if not ok: raise RuntimeError(message)


def floor(x): return x.numerator//x.denominator

def ceil(x): return -floor(-x)

gamma=F(5878847,19531250); D=8*5**10
fee=(F(3,2)+gamma)/5
need(fee==F(17587861,48828125),'uniform fee')
need(3*fee>1,'three-excursion family closes')
Astar=D*(F(3,8)-gamma)
need(Astar.denominator==1 and Astar%8==7,'arithmetic envelope equality clock')
Sstar=(5*Astar-3)/8
need(Sstar.denominator==1,'arithmetic envelope equality integer score')
need(F(5,8)*Astar-F(Astar,5*D)-Sstar==fee,'exact envelope minimum')

patterns=[]
original=[(0,1),(1,1),(2,1),(3,1),(0,2),(1,2),(0,3)]
for a,b in original:
    # 4J+8b gamma < 40-8a; larger J already pays.
    Jmax=ceil(10-2*a-2*b*gamma)-1
    ps=[]
    for js in combinations_with_replacement(range(3,Jmax+1),b):
        J=sum(js)
        if J>Jmax: continue
        bounds=[]
        for j in js:
            hi=floor(D*(F(j,8)-gamma))
            r=5*j%8
            hi-= (hi-r)%8
            first=3+(r-3)%8
            need(first<=hi,'nonempty arithmetic block interval')
            bounds.append({'j':j,'clock_residue_mod8':r,'first_clock':first,'last_clock':hi})
        lower=D*(F(5*J,8)-5+a)
        maxsum=sum(z['last_clock'] for z in bounds)
        need(maxsum>lower,'nonempty total-clock necessary window')
        ps.append({'j_multiset':js,'individual_clock_bands':bounds,
            'sum_nonone_clocks_strict_lower':str(lower),
            'sum_nonone_clocks_upper':maxsum})
    patterns.append({'a':a,'b':b,'defect_partitions':ps})
need([len(z['defect_partitions']) for z in patterns]==[7,5,3,1,4,1,0],'21 partition classification')
# Hostile controls show why mere positive gamma without integrality did not
# close this family, and why the new argument still leaves six open.
need(3*gamma<1,'old margin insufficient')
need(fee+F(3,5)<1 and 2*fee+F(1,5)<1,'six families not declared closed')
print(json.dumps({'status':'PASS','inherited_gamma':str(gamma),'uniform_nonone_fee':str(fee),
 'arithmetic_envelope_equality_only':{'clock':int(Astar),'integer_score':int(Sstar),'actual_word_known':False},
 'closed_family':{'a':0,'b':3,'excess_payment':str(3*fee-1)},
 'remaining_family_count':6,'remaining_defect_partition_count':21,'families':patterns,
 'E12':'OPEN','all_layer_content':'OPEN','universal_Collatz':'OPEN'},indent=2))
