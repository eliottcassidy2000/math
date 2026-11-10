#!/usr/bin/env python3
"""Self-contained exact finite controls; stdlib only. No external data."""
import itertools, json, math
from pathlib import Path

def near(a):
    n=a+2
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    rows[0]^=1<<(n-1); rows[-1]|=1
    return rows

def join(A,B):
    n=len(A); m=len(B)
    return [r|(((1<<m)-1)<<n) for r in A]+[r<<n for r in B]

def hp(A):
    n=len(A)
    if not n:return 1
    dp=[{} for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in dp[mask].items():
            avail=A[v]&~mask
            while avail:
                bit=avail&-avail;avail-=bit;w=bit.bit_length()-1
                d=dp[mask|bit];d[w]=d.get(w,0)+c
    return sum(dp[-1].values())

def hp_perm(A):
    return sum(all(A[u]>>v&1 for u,v in zip(P,P[1:])) for P in itertools.permutations(range(len(A))))

def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def primes(k):
    ans=[];n=2
    while len(ans)<k:
        if all(n%p for p in ans if p*p<=n):ans.append(n)
        n+=1
    return ans

P=primes(12)
def phi(n):
    return math.prod(((1<<(1<<P.index(p)))+1)**e for p,e in factors(n).items())

def encode(n):
    out=[]
    for p,e in factors(n).items():
        for _ in range(e):out=join(out,near(1<<P.index(p)))
    return out

def v2(n):return (n&-n).bit_length()-1

def main():
    rows=[]
    for a in range(1,11):
        x=hp(near(a));assert x==2**a+1
        if a<=6:assert x==hp_perm(near(a))
        rows.append({'a':a,'order':a+2,'H':x})
    tests=0
    for a in range(1,65):
        for b in range(1,65):
            expected=2**math.gcd(a,b)+1 if v2(a)==v2(b) else 1
            assert math.gcd(2**a+1,2**b+1)==expected;tests+=1
    # Actual adjacency verification, not merely assigned scalar values.
    actual=[]
    for n in range(1,17):
        T=encode(n)
        if len(T)<=15:
            h=hp(T);assert h==phi(n)
            actual.append({'n':n,'order':len(T),'H':h})
    # Limit 31 covers first eleven primes without gigantic integers.
    pairs=0
    for a in range(1,32):
        for b in range(1,32):
            A=phi(a);B=phi(b)
            assert phi(a*b)==A*B
            assert (b%a==0)==(B%A==0)
            assert phi(math.gcd(a,b))==math.gcd(A,B)
            assert phi(math.lcm(a,b))==math.lcm(A,B)
            assert (math.gcd(a,b)==1)==(math.gcd(A,B)==1)
            pairs+=1
    optimum=[]
    for k in range(1,6):
        candidates=[c for c in itertools.combinations(range(1,17),k) if len({v2(x) for x in c})==k]
        best=min(sum(x+2 for x in c) for c in candidates)
        winners=[list(c) for c in candidates if sum(x+2 for x in c)==best]
        assert best==2**k-1+2*k
        assert winners==[[2**j for j in range(k)]]
        optimum.append({'k':k,'minimum_total_order':best,'exponent_sets':winners})
    hostile=[join(near(1),near(2)),join(near(1),near(4)),join(near(2),near(4))]
    assert [hp(T) for T in hostile]==[15,51,85]
    values=[15**a*51**b*85**c for a,b,c in itertools.product(range(7),repeat=3)]
    assert len(set(values))==343
    # E_3 embeds into E_4 via vertices 0,1,3, but 3 does not divide 5.
    T=near(2);subset=[0,1,3]
    induced=[sum(1<<j for j,v in enumerate(subset) if T[u]>>v&1) for u in subset]
    assert induced==near(1)
    result={'status':'FINITE-EXACT controls; all-order statements require PROOF.md',
      'near_transitive_H':rows,'gcd_formula_pairs':tests,'arithmetic_pairs_1_to_31':pairs,
      'actual_encoded_tournaments':actual,'bank_optimization_exponents_1_to_16':optimum,
      'hostiles':{'injective_box_0_to_6_weights':[15,51,85],'distinct_values':343,
       'false_divisibility':{'15_divides_51_times_85':True,'quotient':289},
       'contained_E3_in_E4':subset,'addition':{'phi_2':phi(2),'phi_3':phi(3),'phi_5':phi(5)},
       'prime_input_13_composite_image':{'phi_13':phi(13),'factorization':[641,6700417]}}}
    assert 641*6700417==phi(13)
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
