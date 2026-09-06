#!/usr/bin/env python3
"""Exact controls for the conditional-row exponential-rate improvement.

No repository/producer imports. Compare column-mask rook DP with component
matching polynomials, literal column permutations with conditional means,
all small row subsets with the uniform deficit, and positive large examples.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb, factorial, gcd
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def fall(n,k):
    return 0 if k>n else factorial(n)//factorial(n-k)


def parts(n,minimum=2):
    if n==0:
        yield ()
    for a in range(minimum,n+1):
        for tail in parts(n-a,a):
            yield (a,)+tail


def skeleton(partition):
    rows=[]
    offset=0
    for a in partition:
        rows.extend((offset+i,offset+(i+1)%a) for i in range(a))
        offset+=a
    return tuple(tuple(sorted(pair)) for pair in rows)


def reorder(rows,rho):
    out=[None]*len(rows)
    for i,pair in enumerate(rows):
        out[rho[i]]=pair
    return tuple(out)


def defect(rows,sigma):
    counts=Counter(sigma[c]-r for r,pair in enumerate(rows) for c in pair)
    return sum(max(y-2,0) for y in counts.values())


def rook_mask(rows,S):
    counts={0:1}
    for r in S:
        out=counts.copy()
        for mask,count in counts.items():
            for col in rows[r]:
                bit=1<<col
                if not mask&bit:
                    out[mask|bit]=out.get(mask|bit,0)+count
        counts=out
    result=[0]*(len(S)+1)
    for mask,count in counts.items():
        result[mask.bit_count()]+=count
    return result


def multiply(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j]+=x*y
    return out


def rook_components(rows,S):
    n=len(rows)
    adjacency={}
    for r in S:
        for c in rows[r]:
            adjacency.setdefault(r,set()).add(n+c)
            adjacency.setdefault(n+c,set()).add(r)
    visited=set()
    result=[1]
    for root in adjacency:
        if root in visited:
            continue
        visited.add(root)
        todo=[root]
        for x in todo:
            for y in adjacency[x]:
                if y not in visited:
                    visited.add(y)
                    todo.append(y)
        q=sum(len(adjacency[x]) for x in todo)//2
        cyclic=all(len(adjacency[x])==2 for x in todo)
        if cyclic:
            p=[1]+[comb(q-k,k)+comb(q-k-1,k-1) for k in range(1,q//2+1)]
        else:
            need(sum(len(adjacency[x])==1 for x in todo)==2,
                 "noncyclic induced component is a path")
            p=[comb(q+1-k,k) for k in range((q+1)//2+1)]
        result=multiply(result,p)
    need(len(result)==len(S)+1,"induced rook degree equals selected row count")
    return result


def conditional_mean(rows):
    n=len(rows)
    out=Q(0)
    for d in range(1-n,n):
        S=tuple(r for r in range(n) if 0<=r+d<n)
        p=rook_components(rows,S)
        out+=sum(Q((-1)**(k+1)*(k-2)*p[k],fall(n,k))
                 for k in range(3,len(S)+1))
    return out


def conditional_B4(n):
    return Q(2*(n-9),15)+Q(2,n*(n-1))


def main():
    print("ROW-FROZEN DIAGONAL DEFECT: EXACT CONTROLS")
    K=14
    e2_lo=sum(Q(2**k,factorial(k)) for k in range(K+1))
    e2_hi=e2_lo+Q(2**(K+1),factorial(K+1))/(1-Q(2,K+2))
    need(7<e2_lo<e2_hi<Q(37,5),"rational exponential interval")
    alpha_lo=1-5/e2_lo
    alpha_hi=1-5/e2_hi
    need(3*e2_hi-9/e2_hi<21,"conditional lower error constant21")
    need(3*e2_hi+13/e2_hi-4<20,"conditional upper error constant20")

    subset_count=0
    profile_count=0
    for n in range(2,9):
        for partition in parts(n):
            profile_count+=1
            rows=skeleton(partition)
            for mask in range(1<<n):
                S=tuple(i for i in range(n) if mask>>i&1)
                L=len(S)
                p=rook_mask(rows,S)
                need(p==rook_components(rows,S),"mask/component rook equality")
                counts=Counter(c for r in S for c in rows[r])
                double=sum(count==2 for count in counts.values())
                need(double<=L,"selected double-column budget")
                if L>=2:
                    collisions=sum(x==y for r,s in permutations(S,2)
                                   for x in rows[r] for y in rows[s])
                    need(collisions==2*double,"exact ordered collision numerator")
                    need(Q(collisions,4*L*(L-1))<=Q(1,2*(L-1)),
                         "conditional row-pair collision probability")
                theta=Q(L,n)
                for k in range(2,n+4):
                    actual=Q(factorial(k)*p[k],fall(n,k)) if k<=L else Q(0)
                    delta=(2*theta)**k-actual
                    bound=Q(2**k*k*(k-1),n)*theta**(k-1)
                    need(0<=delta<=bound,"uniform conditional all-k deficit")
                if L>=3:
                    exact3=Q(4*fall(L,3),3)-2*double*(L-2)
                    need(p[3]==exact3,"exact induced third matching coefficient")
                subset_count+=1
    print(f"all {profile_count} cycle profiles n2..8, all {subset_count} row subsets PASS")

    literal_count=0
    rho_count=0
    for n in range(2,7):
        P=list(permutations(range(n)))
        for partition in parts(n):
            G=skeleton(partition)
            # Exhaustive row orders through n5; four explicit controls at n6.
            R=P if n<=5 else [P[0],P[1],P[len(P)//3],P[-1]]
            values=[]
            for rho in R:
                rows=reorder(G,rho)
                literal=Q(sum(defect(rows,sigma) for sigma in P),len(P))
                calculated=conditional_mean(rows)
                need(literal==calculated,"literal column law versus conditional rook mean")
                if n>=4:
                    need(calculated>=conditional_B4(n),"conditional fourth-order lower bound")
                    need(n*alpha_hi-21<=calculated<=n*alpha_lo+20,
                         "uniform finite conditional mean bounds")
                values.append(calculated)
                literal_count+=len(P)
                rho_count+=1
            print(f"n={n} type={partition}: row orders={len(R)}, "
                  f"conditional mean range=[{min(values)},{max(values)}]")
            if n==5 and partition==(5,):
                need(Counter(values)==Counter({Q(11,10):40,Q(17,15):40,Q(7,6):40}),
                     "complete first row-dependence hostile")
                need(conditional_mean(reorder(G,(0,1,2,3,4)))==Q(11,10),
                     "hostile identity row labels")
                need(conditional_mean(reorder(G,(0,2,3,1,4)))==Q(7,6),
                     "hostile alternate row labels")
            if n<=4:
                need(len(set(values))==1,"smaller-size row-independence control")
    print(f"literal universe: {rho_count} row orders, {literal_count} column labellings")

    for n,partition,multiplier in [(96,(2,)*48,1),(96,(2,)*48,5),
                                    (97,(97,),1),(97,(97,),5)]:
        need(gcd(n,multiplier)==1,"large-control row labels are a permutation")
        rows=reorder(skeleton(partition),tuple(multiplier*i%n for i in range(n)))
        value=conditional_mean(rows)
        need(n*alpha_hi-21<=value<=n*alpha_lo+20,"large exact conditional envelope")
        need(value>=conditional_B4(n)>0,"positive large fourth-order target")
        need(n*alpha_lo-21>0,"positive asymptotic target envelope")
        print(f"large exact n={n} components={len(partition)} row multiplier={multiplier}: "
              f"floor(mean)={value.numerator//value.denominator}, positive envelope PASS")

    for n in range(4,101):
        raw=Q(2*(n**4-6*n**3+11*n*n-3*n-6),3*fall(n,3))-Q(4*(2*n-3),15)
        need(raw==conditional_B4(n),"summed conditional B4 algebra")
    need(Q(4,225)/8==Q(1,450),"one-permutation asymptotic finite-truncation rate")
    print("one-permutation concentration denominator8(n-1), limiting rate alpha^2/8")
    print(f"PASS {GATES} always-active gates")


if __name__=="__main__":
    main()
