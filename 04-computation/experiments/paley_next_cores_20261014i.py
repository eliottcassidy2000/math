"""Exact injective Paley graph reductions through ten edges, q=3 mod4.

Outputs are polynomial combinations of marked canonical terminal graph sums.
No sampled values are used to infer symbolic formulas.
"""
from dataclasses import dataclass, replace
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial, isqrt, prod
from fractions import Fraction
from collections import Counter
from pathlib import Path
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[2]
POLY = (0, ())


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError("exact integer outside domain")
    return n


def add(a,b):
    out = [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
           for i in range(max(len(a),len(b)))]
    while len(out)>1 and out[-1]==0:
        out.pop()
    return tuple(out)


def scale(a,c):
    return tuple(c*x for x in a)


def linear(a,c):
    return add((0,)+a,scale(a,c))


def merge(result, other, multiplier=1):
    for core,poly in other:
        result[core] = add(result.get(core,(0,)),scale(poly,multiplier))
        if not any(result[core]):
            del result[core]


def normalize(n, edges):
    sign, parity = 1, set()
    for a,b in edges:
        if a==b:
            return 0, ()
        if a>b:
            a,b,sign=b,a,-sign
        edge=(a,b)
        if edge in parity:
            parity.remove(edge)
        else:
            parity.add(edge)
    return sign,tuple(sorted(parity))


@lru_cache(None)
def canonical(n, edges):
    """Exact signed isomorphism key; conflicting signs certify zero."""
    best=None
    signs=set()
    for p in permutations(range(n)):
        sign, target=normalize(n,tuple((p[a],p[b]) for a,b in edges))
        if best is None or target<best:
            best,signs=target,{sign}
        elif target==best:
            signs.add(sign)
    if len(signs)==2:
        return 0,(n,best)
    return next(iter(signs)),(n,best)


def validate(n, edges):
    integer(n)
    if n>20:
        raise ValueError("public receipt scope is at most twenty source vertices")
    if type(edges) is not tuple or len(edges)>10:
        raise ValueError("at most ten exact ordered edge pairs required")
    for edge in edges:
        if type(edge) is not tuple or len(edge)!=2:
            raise ValueError("exact edge pair required")
        for v in edge:
            integer(v)
            if v>=n:
                raise ValueError("vertex outside marked source")
    return n,edges


def evaluate(n, edges):
    validate(n,edges)
    return _evaluate(n,edges)


def _evaluate(n, edges):
    sign,edges=normalize(n,edges)
    if sign==0 or len(edges)%2:
        return ()
    return tuple((core,scale(poly,sign)) for core,poly in reduced(n,edges))


@lru_cache(None)
def reduced(n, edges):
    if n==0:
        return ((POLY,(1,)),)
    neighbors=[[] for _ in range(n)]
    for a,b in edges:
        neighbors[a].append(b)
        neighbors[b].append(a)
    v=min(range(n),key=lambda a:(len(neighbors[a]),-a))
    degree=len(neighbors[v])
    if degree>=3:
        sign,core=canonical(n,edges)
        return ((core,(sign,)),) if sign else ()
    relabel=lambda a:a-(a>v)
    retained=tuple((relabel(a),relabel(b)) for a,b in edges if v not in (a,b))
    if degree==0:
        return tuple((core,linear(poly,-(n-1)))
                     for core,poly in _evaluate(n-1,retained))
    result={}
    if degree==2:
        orient=(-1)**sum(a<v for a in neighbors[v])
        merge(result,_evaluate(n-1,retained),-orient)
    for j in range(n):
        if j==v or j in neighbors[v]:
            continue
        extra=tuple((relabel(j if a==v else a),relabel(j if b==v else b))
                    for a,b in edges if v in (a,b))
        merge(result,_evaluate(n-1,retained+extra),-1)
    return tuple(sorted(result.items()))


def run_source(lengths):
    if type(lengths) is not tuple or any(type(a) is not int or a<1 for a in lengths):
        raise ValueError("positive exact run lengths required")
    if sum(lengths)>10:
        raise ValueError("ten-edge scope")
    n,edges=0,[]
    for length in lengths:
        edges.extend((n+i,n+i+1) for i in range(length))
        n+=length+1
    return n,tuple(edges)


@dataclass(frozen=True)
class Receipt:
    vertices:int
    edges:tuple
    terms:tuple


def compile_receipt(n,edges):
    return Receipt(n,edges,evaluate(n,edges))


def audit_receipt(receipt,n,edges):
    validate(n,edges)
    if type(receipt) is not Receipt or type(receipt.vertices) is not int:
        raise ValueError("typed source receipt required")
    validate(receipt.vertices,receipt.edges)
    if type(receipt.terms) is not tuple:
        raise ValueError("canonical exact term tuple required")
    for term in receipt.terms:
        if type(term) is not tuple or len(term)!=2:
            raise ValueError("canonical term required")
        core,poly=term
        if type(core) is not tuple or len(core)!=2:
            raise ValueError("canonical core required")
        validate(*core)
        if type(poly) is not tuple or any(type(c) is not int for c in poly):
            raise ValueError("integer polynomial required")
    if receipt != compile_receipt(n,edges):
        raise ValueError("wrong source or altered reduction")
    return receipt


def poly_value(poly,q):
    total=0
    for c in reversed(poly):
        total=total*q+c
    return total


def prime(q):
    return type(q) is int and q>=2 and all(q%d for d in range(2,isqrt(q)+1))


def brute(n, edges, q):
    """Even-edge affine-normalized injection census, odd prime q=3mod4."""
    validate(n,edges)
    if not prime(q) or q%4!=3:
        raise ValueError("prime q=3 mod4 required by this finite checker")
    if len(edges)%2:
        return 0
    if n>q:
        return 0
    if n<2:
        return q**n if not edges else 0
    sq={x*x%q for x in range(1,q)}
    chi=[0]+[1 if x in sq else -1 for x in range(1,q)]
    total=0
    for tail in permutations(range(2,q),n-2):
        values=(0,1)+tail
        term=1
        for a,b in edges:
            term*=chi[(values[b]-values[a])%q]
        total+=term
    return q*(q-1)*total


def total_value(terms,q):
    return sum(poly_value(poly,q)*(1 if core==POLY else brute(*core,q))
               for core,poly in terms)


def terminal_census():
    counts={}
    survivors=set()
    for n in range(4,7):
        all_edges=tuple(combinations(range(n),2))
        for m in range((3*n+1)//2,min(10,len(all_edges))+1):
            for edges in combinations(all_edges,m):
                degrees=[0]*n
                for a,b in edges:
                    degrees[a]+=1;degrees[b]+=1
                if min(degrees)<3:
                    continue
                sign,core=canonical(n,edges)
                counts[core]=counts.get(core,0)+1
                if m%2==0 and sign:
                    survivors.add(core)
    return counts,sorted(survivors)


def signed_walk_counts(q,upto):
    """Independent self-avoiding chi-walk subset DP; no graph reduction.

    q<=19 implies every accumulated absolute sum <=q!<2^63.
    The NumPy/Numba arrays use exact signed int64 arithmetic.
    """
    integer(upto)
    if not prime(q) or q%4!=3 or q>19 or upto>=q:
        raise ValueError("finite DP scope: prime q=3mod4<=19 and length<q")
    import numpy as np
    from numba import njit
    sq={x*x%q for x in range(1,q)}
    chi=np.array([0]+[1 if x in sq else -1 for x in range(1,q)],dtype=np.int64)
    @njit
    def kernel(q,upto,chi):
        m=q-1
        size=1<<m
        dp=np.zeros((size,m),dtype=np.int64)
        bits=np.zeros(size,dtype=np.uint8)
        totals=np.zeros(upto+1,dtype=np.int64)
        totals[0]=1
        for v in range(m):
            dp[1<<v,v]=chi[v+1]
        for mask in range(1,size):
            bits[mask]=bits[mask>>1]+(mask&1)
            k=bits[mask]
            if k>upto:
                continue
            for v in range(m):
                value=dp[mask,v]
                totals[k]+=value
                if value and k<upto:
                    for w in range(m):
                        if not (mask&(1<<w)):
                            dp[mask|(1<<w),w]+=value*chi[(w-v)%q]
        return totals*q
    return tuple(map(int,kernel(q,upto,chi)))


def encoded(terms):
    return [{"vertices":core[0],"edges":core[1],"coefficient_ascending":poly}
            for core,poly in terms]


def main():
    checks=0
    def check(condition):
        nonlocal checks
        checks+=1
        if not condition:
            raise ArithmeticError(f"control {checks} failed")
    def reject(call):
        try:
            call()
        except (ValueError,TypeError):
            check(True)
        else:
            raise ArithmeticError("invalid receipt accepted")
    counts,cores=terminal_census()
    check(len(counts)==10 and len(cores)==2)
    wheel=((0,1),(0,2),(0,3),(0,4),(1,2),(2,3),(3,4),(4,1))
    path_complement=tuple((a,b) for a,b in combinations(range(6),2) if b!=a+1)
    sw,ew=normalize(5,wheel)
    cw,wkey=canonical(5,ew)
    cv,vkey=canonical(6,path_complement)
    check(sw*cw==1 and cv==-1)
    check(cores==[wkey,vkey])
    rows=[]
    for lengths in ((2,),(4,),(6,),(8,),(10,)):
        n,edges=run_source(lengths)
        terms=evaluate(n,edges)
        audit_receipt(compile_receipt(n,edges),n,edges)
        check(True)
        rows.append({"runs":lengths,"terms":encoded(terms)})
    expected_A10=(0,-97188,161384,-78252,15334,-1320,42)
    a10=evaluate(*run_source((10,)))
    check(dict(a10)=={POLY:expected_A10,wkey:(-34,)})
    run_partitions=((10,),(8,2),(6,4),(6,2,2),(4,4,2),(4,2,2,2),(2,2,2,2,2))
    total={}
    joints=[]
    for lengths in run_partitions:
        terms=evaluate(*run_source(lengths))
        weight=Fraction(1,prod(factorial(m) for m in Counter(lengths).values()))
        merge(total,terms,weight)
        joints.append({"runs":lengths,"multiset_weight":str(weight),"terms":encoded(terms)})
        for q in (3,7,11):
            check(total_value(terms,q)==brute(*run_source(lengths),q))
    check(total[wkey]==(1166,-173,7) and total[vkey]==(-64,))
    # This direct product uses the stated factored polynomial, not interpolation.
    cubic=(0,3,-4,1) # q(q-1)(q-3)
    degree7=(-32287320,17639622,-4096079,539149,-44330,2308,-71,1)
    explicit=[0]*(len(cubic)+len(degree7)-1)
    for i,a in enumerate(cubic):
        for j,b in enumerate(degree7):
            explicit[i+j]+=a*b
    check(total[POLY]==tuple(Fraction(x,120) for x in explicit))
    # Complete representative cores and source relabellings, independent literal sums.
    core_rows=[]
    for core in counts:
        n,edges=core
        value=evaluate(n,edges)
        for q in (3,7,11):
            check(total_value(value,q)==brute(n,edges,q))
        core_rows.append({"vertices":n,"edges":edges,"labelled_census_count":counts[core],
                          "reduction":encoded(value)})
    check(sum(counts.values())==634)
    # Complete small source universe, not just graphs appearing in the run recursion.
    for n in range(5):
        pairs=tuple(combinations(range(n),2))
        for mask in range(1<<len(pairs)):
            edges=tuple(edge for k,edge in enumerate(pairs) if mask&(1<<k))
            for q in (3,7):
                check(total_value(evaluate(n,edges),q)==brute(n,edges,q))
    # Guard: squaring an edge is1 only on distinct points.
    check(sum(0 if x==y else 1 for x in range(7) for y in range(7))==42!=49)
    check(total_value(evaluate(2,((0,1),(0,1))),7)==42)
    check(total_value(evaluate(2,((0,1),(1,0))),7)==-42)
    check(evaluate(2,((0,0),))==())
    # Retain polynomial factors multiplying a terminal when isolates are deleted.
    wi=evaluate(6,wheel)
    check(dict(wi)=={wkey:(-5,1)})
    check(total_value(wi,7)==brute(6,wheel,7))
    # Raw quotient's Euler parity does not survive degree-two elimination.
    odd_w=sum(sum(v in e for e in wkey[1])%2 for v in range(5))
    odd_v=sum(sum(v in e for e in vkey[1])%2 for v in range(6))
    check(odd_w==odd_v==4)
    check(dict(a10).get(wkey)==(-34,))
    # The global odd-edge involution is a signed cancellation, not no embeddings.
    check(evaluate(3,((0,1),(1,2),(2,0)))==())
    check(7*6*5>0)
    chi5=(0,1,-1,-1,1)
    a2_at5=sum(chi5[(b-a)%5]*chi5[(c-b)%5]
              for a,b,c in permutations(range(5),3))
    check(a2_at5==-20!=total_value(evaluate(*run_source((2,))),7))
    check(poly_value(dict(evaluate(*run_source((2,))))[POLY],5)==20!=-20)
    reject(lambda:brute(*run_source((2,)),5))
    # Reversal of one factor changes the complete value, preserving the source.
    for n,edges in (run_source((4,)),run_source((6,2)),(6,path_complement)):
        reversed_first=((edges[0][1],edges[0][0]),)+edges[1:]
        check(dict(evaluate(n,reversed_first))=={c:scale(p,-1) for c,p in evaluate(n,edges)})
    r=compile_receipt(*run_source((4,)))
    for call in (lambda:evaluate(True,()),lambda:evaluate(21,()),
                 lambda:evaluate(3,((False,1),)),lambda:evaluate(3,((0,3),)),
                 lambda:evaluate(3,[(0,1)]),lambda:run_source((True,)),
                 lambda:run_source((0,)),lambda:run_source((11,)),
                 lambda:audit_receipt(replace(r,vertices=True),*run_source((4,))),
                 lambda:audit_receipt(replace(r,terms=((POLY,(0.,)),)),*run_source((4,))),
                 lambda:audit_receipt(r,*run_source((2,))),
                 lambda:signed_walk_counts(19,19)):
        reject(call)
    # Independent signed-walk DP, including requested q19/length10 test.
    dp_rows=[]
    for q in (7,11,19):
        vals=signed_walk_counts(q,min(q-1,10))
        for k in range(1,len(vals)):
            if k%2:
                check(vals[k]==0)
            else:
                check(vals[k]==total_value(evaluate(*run_source((k,))),q))
        dp_rows.append({"q":q,"A0_through_Amax":vals})
    # Old degree-eight formal run coefficient, independently replayed at q11.
    spec=importlib.util.spec_from_file_location("previous_wheel",ROOT/"04-computation/wheel_run_interface_20261014h.py")
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    hk=old.held_karp_coefficients(11)
    check(total_value(tuple(sorted(total.items())),11)/Fraction(factorial(11),1)==hk[10])
    for lengths in ((2,),(4,),(6,),(8,),(6,2),(4,4),(4,2,2),(2,2,2,2)):
        p,c=old.evaluate(*old.run_edges(lengths))
        own=dict(evaluate(*run_source(lengths)))
        check(own.get(POLY,(0,))==p and own.get(wkey,(0,))==(c,))
    print(json.dumps({"checks":checks,"scope":"all odd prime-power q=3mod4; <=10 edges; public source<=20vertices",
          "terminal_iso_classes":len(counts),"surviving_cores":encoded((c,(1,)) for c in cores),
          "complete_terminal_census":core_rows,"single_runs":rows,"ten_edge_joint_runs":joints,
          "S10_terms":encoded(sorted(total.items())),"named_orientation":{"wheel_edges":wheel,
             "V_path_complement_edges":path_complement,"canonical_V_equals_negative_named_V":True},
          "signed_walk_DP":dp_rows,"R11_degree10":str(hk[10]),
          "A10_q7_cancellation":{"polynomial_part":poly_value(expected_A10,7),
                                  "wheel_U":brute(5,wheel,7),"sum":total_value(a10,7)},
          "reducer_cache":str(reduced.cache_info())},indent=2,default=str))


if __name__=="__main__":
    main()
