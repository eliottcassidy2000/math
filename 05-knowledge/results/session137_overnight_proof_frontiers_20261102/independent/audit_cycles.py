"""Independent class generation, subset/permutation cycles and complement types."""
from tournament_kernel import *
from collections import Counter
import json

def carrier_name(cycles):
    m=len(cycles)
    missing=[(i,j) for i,j in combinations(range(m),2) if not(cycles[i][1]&cycles[j][1])]
    degrees=Counter(v for pair in missing for v in pair)
    ds=sorted(degrees.values())
    if not missing:suffix=''
    elif len(missing)==1:suffix='-e'
    elif len(missing)==2 and ds==[1,1,1,1]:suffix='-2K2'
    elif len(missing)==2 and ds==[1,1,2]:suffix='-P3'
    elif len(missing)==3 and ds==[1,1,2,2]:suffix='-P4'
    else:raise RuntimeError(('unexpected small carrier',m,missing))
    return f'K{m}'+suffix

def run():
    rows=[];names={};class_counts=[];checked=0
    for n,bank in class_bank(7):
        if not n:continue
        class_counts.append(len(bank));strong_count=0;us=[];ms=[]
        for c in sorted(bank):
            a=matrix(n,c)
            if n<3 or not strong(a):continue
            strong_count+=1;cycles=odd_cycles(a);m=len(cycles)
            disjoint=[(i,j) for i,j in combinations(range(m),2) if not(cycles[i][1]&cycles[j][1])]
            require(hp(a)==1+2*m+4*len(disjoint),'independent packing formula through7')
            u=0
            for i,(_,mask) in enumerate(cycles):
                universal=all(i==j or (mask&other) for j,(_,other) in enumerate(cycles))
                rest=[j for j in range(n) if not(mask>>j&1)]
                transitive=not any(a[x][y] and a[y][z] and a[z][x]
                                   for x,y,z in permutations(rest,3))
                require(universal==transitive,'universal complement iff');u+=universal
            require(m>=((n-1)**2)//4,'Moon floor')
            require(u>=(1 if n==3 else 2 if n%2==0 else 4),'universal floor')
            ms.append(m);us.append(u);checked+=1
            if m<=11:
                name=carrier_name(cycles)
                names.setdefault(name,{'n':n,'code':c,'H':hp(a),'m':m,
                                       'cycles':[list(p) for p,_ in cycles]})
        rows.append({'n':n,'classes':len(bank),'strong_classes':strong_count,
                     'min_odd_cycles':min(ms,default=None),'min_universal':min(us,default=None)})
    require(class_counts==[1,1,2,4,12,56,456],'class counts')
    expected={'K1','K2','K4','K5','K6','K6-e','K7','K8','K9','K9-e','K9-P4','K10-2K2','K10-P3','K11','K11-e'}
    require(set(names)==expected,'complete small carrier list')
    print(json.dumps({'scope':'FINITE-EXACT all tournament classes through7; complete connected carrier m<=11 conditional only on cited Moon bound and inherited SCC theorem',
                      'summary':rows,'strong_classes_audited':checked,'connected_carriers':names,
                      'hostile_missing_carriers':['K3','P3','K4-e','K8-e']},sort_keys=True,indent=2))
if __name__=='__main__':run()
