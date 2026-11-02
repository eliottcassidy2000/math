"""Independent complete cycle-budget 15 census via all classes through order8.

Cycle enumeration is by subsets/permutations, not the producer's directed DFS.
Only complement edges and degree-preserving permutations are used for graph
isomorphism. Moon's floor excludes source order>=9 when cycle count<=15.
"""
from tournament_kernel import *
from collections import Counter
from math import isqrt
import json

def graph_key(cycles):
    m=len(cycles)
    b=[[int(i!=j and not(cycles[i][1]&cycles[j][1])) for j in range(m)] for i in range(m)]
    vs=[i for i in range(m) if any(b[i])]
    return (m-len(vs),len(vs),canon(induced(b,vs)) if vs else 0)

def run():
    bank={};summary=[]
    for n,reps in class_bank(8):
        if not n:continue
        strong_count=0;qualifying=0
        for code in sorted(reps):
            a=matrix(n,code)
            if n<3 or not strong(a):continue
            strong_count+=1;cycles=odd_cycles(a);m=len(cycles)
            require(m>=((n-1)**2)//4,'Moon floor')
            if m>15:continue
            qualifying+=1;key=graph_key(cycles)
            a2=sum(not(cycles[i][1]&cycles[j][1]) for i,j in combinations(range(m),2))
            h=hp(a);require(h==1+2*m+4*a2,'packing formula')
            bank.setdefault(key,{'key':key,'m':m,'alpha2':a2,'H':h,'n':n,'code':code})
        summary.append({'n':n,'classes':len(reps),'strong_classes':strong_count,'budget15_qualifying_classes':qualifying})
    require([r['classes'] for r in summary]==[1,1,2,4,12,56,456,6880],'class counts')
    expected={(1,0,0),(2,0,0),(4,0,0),(5,0,0),(4,2,1),(6,0,0),(7,0,0),(8,0,0),
              (5,4,44),(7,2,1),(9,0,0),(6,4,12),(7,3,6),(9,2,1),(11,0,0),
              (4,8,189861952),(8,4,12),(10,2,1),(12,0,0),(10,3,6),(11,2,1),(13,0,0),
              (6,8,297024),(8,6,15632),(10,4,12),(10,4,44),(12,2,1),(9,6,22800),(15,0,0)}
    require(set(bank)==expected,'29 carrier types')
    pair_min=[0]+[10**9]*15
    for m in range(1,16):
        for row in bank.values():
            k=row['m']
            if k<=m:pair_min[m]=min(pair_min[m],row['alpha2']+pair_min[m-k]+k*(m-k))
    require(pair_min[3]==2 and pair_min[10]==2 and pair_min[14]==1,'link pair minima')
    minimum31=min(row['m'] for row in bank.values() if row['H']==31)
    require(minimum31==9,'H31 minimal cycle budget')
    print(json.dumps({'scope':'FINITE-EXACT complete all source classes through8; complete connected cycle budget<=15 using cited Moon bound',
                      'summary':summary,'connected_types':sorted(bank.values(),key=lambda r:(r['m'],r['key'])),
                      'minimum_disjoint_pairs_by_cycle_budget':pair_min,'H31_minimal_cycle_budget':minimum31},sort_keys=True,indent=2))
if __name__=='__main__':run()
