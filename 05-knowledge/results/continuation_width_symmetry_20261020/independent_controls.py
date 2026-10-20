#!/usr/bin/env python3
"""Separate exact representations; this is not an independent human/agent audit."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import comb, factorial, gcd
from pathlib import Path
import hashlib
import json
import verify as v

CHECKS = Counter()
smooth = lru_cache(None)(v.smooth_part)


def need(value, label):
    CHECKS[label] += 1
    if not value:
        raise RuntimeError(label)


def kernel_dp(n):
    """Store smooth-kernel multisets rather than future gcd-response vectors."""
    lengths = sorted(range(1,n+1,2),key=lambda a:(v.largest_prime(a),a),reverse=True)
    states = {(0,()):Q(1)}
    peak=1
    for pos,ell in enumerate(lengths):
        next_prime = v.largest_prime(lengths[pos+1]) if pos+1<len(lengths) else 1
        new=defaultdict(Q)
        for (degree,kernels),weight in states.items():
            interaction=sum(gcd(a,ell) for a in kernels)
            for k in range((n-degree)//ell+1):
                exponent=k*(ell-1)//2+ell*comb(k,2)+k*interaction
                factor=Q(2**exponent,ell**k*factorial(k))
                merged=tuple(sorted(smooth(a,next_prime)
                                    for a in kernels+(ell,)*k))
                new[degree+k*ell,merged]+=weight*factor
        states=dict(new)
        peak=max(peak,len(states))
    coefficients=[Q(0)]*(n+1)
    for (d,kernels),weight in states.items():
        coefficients[d]+=weight
    return coefficients,peak


def mul(a,b,n):
    out=[0]*(n+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:n+1-i]):
            out[i+j]+=x*y
    return out


def cycle_type(perm, domain):
    unseen=set(domain)
    sizes=[]
    while unseen:
        start=min(unseen)
        x=start
        size=0
        while x in unseen:
            unseen.remove(x)
            size+=1
            x=perm[x]
        need(x==start,"permutation cycles remain in marked domain")
        sizes.append(size)
    return tuple(sorted(sizes))


def check_group_factorization(adj,order,group,delta):
    n=len(adj)
    orbits=sorted({tuple(sorted({p[v] for p in group})) for v in range(n)})
    owner={v:i for i,o in enumerate(orbits) for v in o}
    pos={v:i for i,v in enumerate(order)}
    marked_orbits={i for i,o in enumerate(orbits) if len(o)>=5}
    for a in range(n):
        for b in range(n):
            if adj[a]>>b&1 and pos[a]>pos[b] and owner[a]!=owner[b]:
                marked_orbits.update((owner[a],owner[b]))
    B={v for i in marked_orbits for v in orbits[i]}
    quiet=[o for i,o in enumerate(orbits) if i not in marked_orbits and len(o)==3]
    t=len(quiet)
    u=n-len(B)-3*t
    H=[p for p in group if all(p[v]==v for v in range(n) if v not in B)]
    need(len(B)<=2*delta,"factorization uses bounded skeleton")
    need(len(group)==len(H)*3**t,"full symmetry product order")
    predicted=Counter()
    for perm in H:
        base=cycle_type(perm,B)
        for j in range(t+1):
            typ=tuple(sorted(base+(3,)*j+(1,)*(u+3*(t-j))))
            predicted[typ]+=comb(t,j)*2**j
    actual=Counter(cycle_type(p,range(n)) for p in group)
    need(actual==predicted,"full permutation cycle index factors")
    for ports in [()] + [(i,) for i in range(n)] + list(combinations(range(n),2)):
        h=sum(bool(set(o)&set(ports)) for o in quiet)
        predicted_order=sum(all(p[v]==v for v in ports if v in B)
                            for p in H)*3**(t-h)
        actual_order=sum(all(p[v]==v for v in ports) for p in group)
        need(predicted_order==actual_order,"one and two marked port stabilizers factor")


def atom_words(n):
    if n==0:
        yield ()
    for a in (1,3):
        if n>=a:
            for rest in atom_words(n-a):
                yield (a,)+rest


def ordered_join(blocks):
    n=sum(map(len,blocks))
    out=[0]*n
    start=0
    for block in blocks:
        after=sum(1<<j for j in range(start+len(block),n))
        for i,row in enumerate(block):
            out[start+i]=(row<<start)|after
        start+=len(block)
    return out


def atom_word(word):
    return ordered_join([v.adjacency(1,0) if a==1 else
                         v.adjacency(3,1<<1) for a in word])


def low_slack_samples():
    samples=0
    for n in range(2,8):
        for word in atom_words(n):
            U=atom_word(word)
            central=ordered_join([[0],U,[0]])
            central[0]^=1<<(len(central)-1)
            central[-1]^=1
            for left,right in [((),()),((1,),(3,)),((3,),(1,))]:
                adj=ordered_join([atom_word(left),central,atom_word(right)])
                f,s,d,_,aut=v.check_tournament(adj)
                triangles=(word+left+right).count(3)
                need((f,s,d,aut)==(triangles+1,3*triangles,3,3**triangles),
                     "all generated first-slack words have exact invariants")
                samples+=1
    # Explicit diagonal triangle synchronization hostile.
    adj=[0]*6
    for i in range(3):
        adj[i]|=1<<((i+1)%3)
        adj[3+i]|=1<<(3+(i+1)%3)
        for j in range(3):
            if i==j: adj[i]|=1<<(3+j)
            else: adj[3+j]|=1<<i
    f,order=v.feedback(adj)
    group=v.automorphisms(adj)
    need(f==4 and len(group)==3,"synchronized triangle exact feedback and group")
    need(all(p[3+i]==3+p[i] for p in group for i in range(3)),
         "triangle phases are locked")
    need(all({p[i] for p in group}==set(range(3)) for i in range(3)),
         "first synchronized orbit transitive")
    check_group_factorization(adj,order,group,6)
    return {"first_slack_samples":samples,
            "synchronized_hostile":{"adjacency":adj,"feedback":f,"order":order,
                                    "automorphisms":group}}


def main():
    kernels=[]
    for n in (12,24,40,60):
        coeff,peak=kernel_dp(n)
        expected,info=v.response_dp(n,"prime")
        need(coeff==expected,"smooth-multiset DP equals response DP")
        kernels.append({"N":n,"smooth_multiset_peak":peak,"response_peak":info["peak"]})
    # Explicit subset-zeta Gram factorization, without LDL elimination.
    ps=(7,11,13,17)
    masks=range(1<<len(ps))
    for a in masks:
        for b in masks:
            gram=sum(
                __import__("functools").reduce(
                    lambda x,i:x*(2**(ps[i]-1)-1),
                    (i for i in range(len(ps)) if c>>i&1),1)
                for c in masks if c&a==c and c&b==c)
            direct=2**sum(ps[i]-1 for i in range(len(ps)) if (a&b)>>i&1)
            need(gram==direct,"subset-zeta factorization exact")
    # Verify symmetry order and class bounds from all labelled small objects.
    class_rows={}
    for n in range(7):
        by_delta=defaultdict(Q)
        group_sizes=Counter()
        for mask in range(1<<comb(n,2)):
            adj=v.adjacency(n,mask)
            f,order=v.feedback(adj)
            group=v.automorphisms(adj)
            s=sum(any(p[i]!=i for p in group) for i in range(n))
            delta=3*f-s
            need(len(group)<=3**f,"group order paid by feedback")
            need((len(group)==3**f)==(delta==0),"group bound equality exact")
            check_group_factorization(adj,order,group,delta)
            by_delta[delta]+=Q(len(group),factorial(n))
            group_sizes[len(group)]+=1
        need(all(x.denominator==1 for x in by_delta.values()),"class census integral")
        base=[0]*(n+1)
        for j in range(n+1):
            base[j]=(j==0)+(base[j-1] if j else 0)+(base[j-3] if j>=3 else 0)
        square=mul(base,base,n)
        central=[0]*(n+1)
        for j in range(2,n+1):
            central[j]=base[j-2]-(j==2)-(j==3)
        exact3=mul(central,square,n)[n]
        exact4=square[n-5] if n>=5 else 0
        need(by_delta.get(3,0)==exact3,"first nonzero slack exact class series")
        need(by_delta.get(4,0)==exact4,"five-orbit slack exact class series")
        need(all(not by_delta.get(d,0) for d in (1,2,5)),"low absent slack values")
        upper=[0]*(n+1)
        power=[1]+[0]*n
        bounds={}
        for m in range(9):
            power=mul(power,base,n)
            for j in range(m,n+1):
                upper[j]+=2**comb(m,2)*power[j-m]
            if m%2==0:
                D=m//2
                actual=sum(x for delta,x in by_delta.items() if delta<=D)
                need(base[n]<=actual<=upper[n],"fixed slack generating upper and lower bounds")
                bounds[D]={"actual_classes":int(actual),"upper":upper[n]}
        class_rows[n]={"by_delta":{d:int(c) for d,c in sorted(by_delta.items())},
                       "group_sizes":dict(sorted(group_sizes.items())),
                       "bounds":bounds}
    # Exact finite controls for the elementary odd-factorial induction.
    fact=1
    for n in range(1,202):
        fact*=n
        if n>=3 and n%2:
            odd=fact
            while odd%2==0:
                odd//=2
            need(odd<=3**((n*n-1)//8),"odd factorial entropy bound")
            need((odd==3**((n*n-1)//8))==(n==3),"odd factorial strict boundary")
            if n>=5:
                need((n+1)*(n+2)//2<=3**((n+1)//2),"inductive multiplier inequality")
    # Port counterexample: clone core alone misses internal triangle position.
    triangle=v.adjacency(3,1<<1)  # 0->1->2->0
    group=v.automorphisms(triangle)
    need(len(group)==3,"triangle rotations")
    need(any(p[0]!=0 for p in group),"unmarked triangle erases selected port")
    need(sum(p[0]==0 for p in group)==1,"one marked triangle port fixes phase")
    low_slack=low_slack_samples()
    print(json.dumps({"status":"FINITE-EXACT, separate representations by same agent",
                      "checks":dict(CHECKS),"total_checks":sum(CHECKS.values()),
                      "kernel_dp":kernels,"class_census":class_rows,"low_slack":low_slack,
                      "main_checker_sha256":hashlib.sha256(
                          Path(v.__file__).read_bytes()).hexdigest()},
                     indent=2,sort_keys=True))


if __name__=="__main__":
    main()
