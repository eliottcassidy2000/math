"""Exact private-certificate incompatibility in the two-step pairing family.
Reproduce with python3 (or python3 -O) on this file. No external packages.
"""
from fractions import Fraction
from itertools import product
import json


def require(condition,message):
    if not condition: raise AssertionError(message)


def pair(n): return (n+1)//2


def step(n,bits):
    i=pair(n)
    flipped=bits.get(i,0)
    if n&1: return i-1 if flipped else 3*i-1
    return 3*i if flipped else i


def private_plans(n):
    answers=[]
    def visit(v,bits,path):
        if len(path)>1 and v<n:
            cost=sum((Fraction(n,i) for i,b in bits.items() if b),Fraction(0))
            answers.append((cost,tuple(sorted(bits.items())),tuple(path)))
            return
        if len(path)==3:return
        i=pair(v)
        for bit in [0,1]:
            if i in bits and bits[i]!=bit:continue
            new=bits|{i:bit}
            y=step(v,new)
            visit(y,new,path+[y])
    visit(n,{},[n])
    return answers


def valuation3(n):
    result=0
    while n%3==0:
        n//=3;result+=1
    return result


def independent_path_cover(limit):
    vertices=set(range(2,limit+1,2))
    edges={(4*k,6*k) for k in range(1,limit//6+1)}
    selected={i for i in vertices if valuation3(i)%2==1}
    matching={(4*k,6*k) for k in range(1,limit//6+1) if valuation3(k)%2==0}
    require(all(a in selected or b in selected for a,b in edges),"proposed set fails cover")
    flattened=[v for edge in matching for v in edge]
    require(len(set(flattened))==len(flattened),"matching edges meet")
    require(len(matching)==len(selected),"cover and matching sizes differ")
    # Independent exhaustive search on each tiny connected path.
    adjacency={v:set() for v in vertices}
    for a,b in edges:adjacency[a].add(b);adjacency[b].add(a)
    unseen=set(vertices)
    optimal=0
    while unseen:
        start=min(unseen)
        component={start};stack=[start]
        while stack:
            for v in adjacency[stack.pop()]:
                if v not in component:component.add(v);stack.append(v)
        unseen-=component
        ordered=sorted(component)
        local_edges=[(a,b) for a,b in edges if a in component]
        best=len(ordered)
        for choice in product([0,1],repeat=len(ordered)):
            take={v for v,b in zip(ordered,choice) if b}
            if all(a in take or b in take for a,b in local_edges):best=min(best,len(take))
        optimal+=best
    require(optimal==len(selected),"independent optimum mismatch")
    return len(vertices),len(selected)


def main():
    private_checks=0
    for n in range(3,4001,4):
        plans=private_plans(n)
        best=min(cost for cost,_,_ in plans)
        minimizers=[p for p in plans if p[0]==best]
        i=pair(n);j=3*i//2
        require(len(minimizers)==1,"private optimum not unique")
        require(minimizers[0][1]==((i,0),(j,1)),"wrong private optimum")
        require(best==Fraction(n,j),"wrong private cost")
        require(minimizers[0][2]==(n,3*i-1,j-1),"wrong private path")
        private_checks+=1
    print(json.dumps({"private_optima_checked":private_checks,"first_source":3,"last_source":3999}))
    for r in range(1,1001):
        n,m=8*r-1,12*r-1
        cn=min(private_plans(n));cm=min(private_plans(m))
        bn,bm=dict(cn[1]),dict(cm[1])
        require(bn[6*r]==1 and bm[6*r]==0,"expected incompatibility absent")
    print(json.dumps({"incompatible_optimum_families_checked":1000,"first_pair":[7,11],"conflict_pair_index":6}))
    for limit in range(2,301):independent_path_cover(limit)
    print(json.dumps({"exact_cover_matching_audits":299,"max_pair_index":300}))
    for limit in [112,223,1000,10000,1000000]:
        if limit<=300:count,cover=independent_path_cover(limit)
        else:
            count=limit//2
            cover=sum((limit//(6*3**(2*j))-limit//(6*3**(2*j+1))) for j in range(20))
        print(json.dumps({"pair_index_limit":limit,"bad_source_policies":count,"minimum_abandoned_optima":cover,
                          "abandoned_fraction":str(Fraction(cover,count))}))
    seed=223
    chain=[]
    n=seed
    while n%4==3:
        optimum=min(private_plans(n));i=pair(n)
        chain.append({"source":n,"pair":i,"v3_pair":valuation3(i),"private_flip":3*i//2,
                      "private_path":optimum[2],"cost":str(optimum[0])})
        n=(3*n+1)//2
    print(json.dumps({"seed":seed,"private_conflict_chain":chain,"first_source_outside_bad_class":n},sort_keys=True))
    # Same saved bits support both individually compatible policies.
    bits={i:valuation3(i)%2 for i in range(1,10000)}
    for n in range(3,10000,4):
        i=pair(n)
        if valuation3(i)%2==0:
            j=3*i//2
            require(bits[i]==0 and bits[j]==1,"compatible policy family not realized")
            require(step(step(n,bits),bits)<n,"selected policy fails")
        else:
            require(step(n,bits)<n,"abandoned optimum was not rescued by own flip")
    for limit in [112,223,1000]:
        private=Fraction(0);common=Fraction(0);gap=Fraction(0)
        for i in range(2,limit+1,2):
            n=2*i-1;j=3*i//2
            private+=Fraction(n,j)
            trace=Fraction(n,i) if bits[i] else Fraction(n,j)
            common+=trace
            if bits[i]:gap+=Fraction(n,3*i)
        require(common-private==gap,"trace gap identity failed")
        print(json.dumps({"trace_pair_cutoff":limit,"source_normalizer":2*limit,
                          "private_mean":float(private/(2*limit)),"common_trace_mean":float(common/(2*limit)),
                          "asymptotic_private_mean":"1/3","asymptotic_common_trace_mean":"3/8"}))
    # This common assignment is not asserted to be globally P_2.
    require(step(step(6,bits),bits)>=6,"control should expose missing good-source constraints")
    print(json.dumps({"compatible_family_assignment":"epsilon_i=v3(i) mod2","global_P2_hostile_source":6,
                      "hostile_path":[6,step(6,bits),step(step(6,bits),bits)]}))
    print("ALL CHECKS PASSED")


if __name__=="__main__":main()
