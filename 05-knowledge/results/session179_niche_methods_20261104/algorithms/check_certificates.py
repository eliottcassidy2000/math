#!/usr/bin/env python3
"""Exact standalone checks; no dependencies; assertions remain live under -O."""
import itertools, json
from pathlib import Path

def check(p, msg):
    if not p: raise RuntimeError(msg)

def cat(u,v):
    A,B,C=u; D,E,F=v
    return (D*A,E*B,D*C+B*F)

def power(u,k):
    out=(1,1,0)
    while k:
        if k&1: out=cat(out,u)
        u=cat(u,u); k//=2
    return out

def compile_word(w,multiplier=3):
    out=(1,1,0)
    for b in w: out=cat(out,(multiplier if b else 1,2,b))
    return out

def compile_runs(w):
    out=(1,1,0)
    for b,group in itertools.groupby(w):
        out=cat(out,power((3 if b else 1,2,b),sum(1 for _ in group)))
    return out

def first_certified(t):
    A,B,C=t
    if A>=B: return None
    r=(-C*pow(A,-1,B))%B
    k=max(0,(C-(B-A)*r)//((B-A)*B)+1)
    return r+B*k

def certified(t,n):
    """Requires t authenticated by the word compiler; arbitrary triples are not certificates."""
    A,B,C=t
    return n>0 and (A*n+C)%B==0 and (B-A)*n>C

def orbit(n,L,multiplier=3):
    bits=[]
    for _ in range(L):
        b=n%2; bits.append(b); n=(multiplier*n+1)//2 if b else n//2
    return bits,n

def main():
    words=0; lift_checks=0; first_checks=0
    for L in range(1,14):
        for word in itertools.product((0,1),repeat=L):
            t=compile_word(word); A,B,C=t
            check(t==compile_runs(word),'run compiler')
            r=(-C*pow(A,-1,B))%B
            for k in (0,1,3):
                n=r+B*k; bits,y=orbit(n,L)
                check(bits==list(word),'residue legality')
                check(B*y==A*n+C,'affine endpoint')
                check(certified(t,n)==(n>0 and y<n),'certificate equivalence')
                lift_checks+=1
            n=first_certified(t)
            if n is not None:
                check(certified(t,n),'first witness acceptance')
                check(n-B<=0 or not certified(t,n-B),'first witness minimality')
                check(orbit(n,L)[1]<n,'first witness actual descent')
                first_checks+=1
            words+=1
    records=[]; max_depth=(0,0)
    for start in range(2,10001):
        x=start; w=[]
        while True:
            b=x%2; w.append(b); x=(3*x+1)//2 if b else x//2
            if x<start: break
        t=compile_runs(w)
        check(certified(t,start),'finite first-descent census')
        if len(w)>max_depth[1]: max_depth=(start,len(w)); records.append(max_depth)
    conservative_max=(0,0)
    for start in range(2,10001):
        x=start; ones=0; w=[]
        while True:
            b=x%2; w.append(b); ones+=b
            x=(3*x+1)//2 if b else x//2
            if x<start and 5*ones<=3*len(w): break
        t=compile_runs(w)
        check(certified(t,start),'conservative grammar certificate')
        A,B,C=t
        check(A**5 * 32**len(w) <= B**5 * 27**len(w),'conservative uniform envelope')
        if len(w)>conservative_max[1]: conservative_max=(start,len(w))
    near=[]
    for a in range(1,501):
        A=3**a; L=A.bit_length(); B=2**L
        check(A<B<2*A,'near-critical exact bound')
        t=cat(power((3,2,1),a),power((1,2,0),L-a))
        check(t[:2]==(A,B),'near-critical compiler')
        n=first_certified(t)
        check(certified(t,n),'near-critical positive certificate')
        if a in (1,2,10,100,500):
            near.append({'ones':a,'length':L,'first_witness_bits':n.bit_length(),'A_bits':A.bit_length(),'C_bits':t[2].bit_length()})
    # Real positive cycles are guards against dropping the strict threshold.
    t=compile_word([1,0]); check(t==(3,4,1),'1 cycle')
    check(not certified(t,1) and first_certified(t)==5,'strict threshold at 1')
    w,y=orbit(13,7,5); t=compile_word(w,5)
    check(y==13 and t==(125,128,39),'5n+1 cycle')
    check(not certified(t,13),'reject multiplier-5 cycle equality')
    # Sanity-check a semilinear family whose base/period slopes all contract.
    # Parikh coordinates (ones, zeros); base (1,1), periods (1,1),(1,2),(0,1).
    semilinear=0
    for x,y,z in itertools.product(range(16),repeat=3):
        a=1+x+y; L=2+2*x+3*y+z
        check(3**a<2**L,'semilinear contraction')
        check(3**(2*a)<=3**L,'uniform sqrt(3)/2 envelope')
        semilinear+=1
    result={'status':'PASS','all_binary_words_lengths_1_to_13':words,
        'independent_positive_and_zero_lift_checks':lift_checks,'minimal_positive_descent_witnesses':first_checks,
        'first_descent_starts_2_to_10000':9999,'first_descent_records':records,
        'conservative_grammar_starts_2_to_10000':9999,'conservative_max_depth':conservative_max,
        'near_critical_exact_checks':500,'near_critical_examples':near,
        'semilinear_uniform_margin_checks':semilinear,
        'hostile_controls':['all-odd expanding words reject','1 cycle equality rejects','5n+1 13 cycle equality rejects'],
        'scope':'Finite checks validate identities and implementation, not Parikh theorem or universal Collatz coverage.'}
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
