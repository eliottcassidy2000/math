"""Exact once-per-pair P2 cut, tree optimum, and a global P4 hostile map.
Run with python3 or python3 -O. Standard library only.
"""
from functools import lru_cache
from itertools import product
from math import gcd
import json


def require(test,message):
    if not test:raise AssertionError(message)


def v3(n):
    answer=0
    while n%3==0:n//=3;answer+=1
    return answer


def clauses(X):
    equal=[];less=[]
    for i in range(2,X+1):
        if i%2==0:
            j=3*i//2
            if j<=X:equal.append((i,j))
        else:
            j,k=(3*i-1)//2,(3*i+1)//2
            if j<=X:less.append((j,i))
            if k<=X:less.append((i,k))
    return equal,less


def feasible(bits,X):
    eq,le=clauses(X)
    return all(bits[i]+bits[j]==1 for i,j in eq) and all(bits[i]<=bits[j] for i,j in le)


def tree_optimum(X):
    zero=[0]*(X+1);one=[0]*(X+1)
    for i in range(X,1,-1):
        a,b=0,1
        if i%2==0:
            j=3*i//2
            if j<=X:a+=one[j];b+=zero[j]
        else:
            j,k=(3*i-1)//2,(3*i+1)//2
            if j<=X:a+=zero[j];b+=min(zero[j],one[j])
            if k<=X:a+=min(zero[k],one[k]);b+=one[k]
        zero[i]=a;one[i]=b
    bits=[0]*(X+1)
    if X>=2:
        stack=[(2,int(one[2]<zero[2]))]
        while stack:
            i,b=stack.pop();bits[i]=b
            if i%2==0:
                j=3*i//2
                if j<=X:stack.append((j,1-b))
            else:
                j,k=(3*i-1)//2,(3*i+1)//2
                if j<=X:stack.append((j,0 if not b else int(one[j]<zero[j])))
                if k<=X:stack.append((k,1 if b else int(one[k]<zero[k])))
    require(feasible(bits,X),'reconstructed optimum violates clause')
    answer=min(zero[2],one[2]) if X>=2 else 0
    require(sum(bits)==answer,'reconstructed objective mismatch')
    return answer,bits


def baseline(X):
    result=0;power=1
    while 3*power<=X:
        result+=X//(3*power)-X//(9*power);power*=9
    return result


def extra(X):return max(0,(X-23)//54+1)


@lru_cache(None)
def original_bit(i):
    if i<3:return 0
    k,r=divmod(i,3)
    return 1-original_bit(2*k) if r==0 else 0 if r==1 else original_bit(2*k+1)


def step(n,edited=False):
    i=(n+1)//2;b=0 if edited and i==23 else original_bit(i)
    return n+i if (n%2)^b else n-i


def first_descent(n,edited=False,horizon=4):
    path=[n]
    for _ in range(horizon):
        path.append(step(path[-1],edited))
        if path[-1]<n:return path
    return None


def main():
    brute_cases=0
    for X in range(2,17):
        best=X
        for choice in product([0,1],repeat=X-1):
            bits=[0,0]+list(choice)
            if feasible(bits,X):best=min(best,sum(choice))
            brute_cases+=1
        require(tree_optimum(X)[0]==best,'tree/brute discrepancy')
    print(json.dumps({'independent_exhaustive_assignments':brute_cases,'cutoffs_checked':[2,16]}))
    first=None
    for X in range(2,224):
        opt,_=tree_optimum(X)
        require(opt>=baseline(X)+extra(X),'new cut exceeds exact optimum')
        if first is None and opt>baseline(X):first=X
    require(first==23,'first strict prefix improvement changed')
    print(json.dumps({'first_prefix_above_matching':first,'matching':baseline(first),'exact_optimum':tree_optimum(first)[0]}))
    for X in [23,33,112,223,1000,3000,10000,100000]:
        opt,_=tree_optimum(X)
        print(json.dumps({'pair_cutoff':X,'exact_P2_prefix_optimum':opt,'new_explicit_lower_bound':baseline(X)+extra(X),
                          'optimum_fraction':opt/X}))
    # Local implication proof checked independently over all six Boolean bits.
    local=0
    for i,j,h,l,a,b in product([0,1],repeat=6):
        if j+l==1 and j<=i<=h<=a and l<=b:
            require(a+b>=1,'local extra cut fails');local+=1
    for s in range(1000):
        i,j,h,l,a,b=16*s+7,24*s+10,24*s+11,36*s+15,36*s+17,54*s+23
        require(j==(3*i-1)//2 and h==(3*i+1)//2 and l==3*j//2,'middle indices')
        require(a==(3*h+1)//2 and b==(3*l+1)//2,'terminal indices')
        require(all(x%2==1 for x in [i,h,l,a,b]) and j%2==0,'parity lost')
        require(v3(a)==v3(b)==0,'terminal is not in unmatched bank')
    require((23-17)%gcd(36,54)!=0,'extra arithmetic progressions intersect')
    X=100000
    matching=[(2*r,3*r) for r in range(1,X//3+1) if v3(r)%2==0]
    endpoints=[x for edge in matching for x in edge]
    added=[(36*s+17,54*s+23) for s in range(extra(X))]
    more=[x for edge in added for x in edge]
    require(len(set(endpoints))==len(endpoints),'canonical matching not disjoint')
    require(len(set(more))==len(more) and set(endpoints).isdisjoint(more),'extra edge collision')
    require(len(matching)==baseline(X),'matching floor count')
    print(json.dumps({'feasible_six_bit_assignments':local,'affine_packets_checked':1000,
                      'matching_edges_at_100000':len(matching),'disjoint_extra_edges_at_100000':len(added),
                      'proved_lower_density':'29/108'}))
    require(original_bit(17)==0 and original_bit(23)==1,'wrong surgery bits')
    delayed=[]
    for n in range(3,47):
        require(first_descent(n,False,2) is not None,'original member fails P2')
        path=first_descent(n,True,4)
        require(path is not None,'edited member fails finite P4 audit')
        if len(path)>3:delayed.append(path)
    require(delayed==[[30,45,68,34,17]],'wrong delayed-source set')
    print(json.dumps({'global_P4_surgery':'unflip only pair23 in free-zero P2 member',
                      'finite_sources_needed_for_global_audit':[3,46],'only_delayed_path':delayed[0],
                      'extra_cut_hostile_bits':{'17':0,'23':0}}))
    print('ALL CHECKS PASSED')


if __name__=='__main__':main()
