#!/usr/bin/env python3
"""Independent exact checks. Does not import either producer or its outputs."""
import itertools, json, math
from functools import lru_cache
from pathlib import Path

def near(a):
    return tuple(tuple((u<v) ^ ({u,v}=={0,a+1}) if u!=v else False for v in range(a+2)) for u in range(a+2))

def fan(n):
    if n==1:return near(1)
    k=2*n; losses={n,n+1,2*n+1}
    return tuple(tuple(False if u==v else (v not in losses if u==0 else (u in losses if v==0 else u<v)) for v in range(k+2)) for u in range(k+2))

def count(A):
    N=len(A)
    @lru_cache(None)
    def tails(last,unvisited):
        if not unvisited:return 1
        return sum(tails(v,unvisited^(1<<v)) for v in range(N) if (unvisited>>v)&1 and A[last][v])
    return sum(tails(v,((1<<N)-1)^(1<<v)) for v in range(N))

def factors(n):
    d={};p=2
    while p*p<=n:
        while n%p==0:d[p]=d.get(p,0)+1;n//=p
        p+=1
    if n>1:d[n]=1
    return d

def v2(n):return (n&-n).bit_length()-1

def main():
    near_rows=[]
    for a in range(1,10):
        h=count(near(a));assert h==2**a+1
        near_rows.append([a,a+2,h])
    fans=[]
    for n in range(1,7):
        A=fan(n);h=count(A);assert h==4**n-1
        # Displayed spanning directed cycle certifies strong connectivity.
        assert all(A[u][(u+1)%len(A)] for u in range(len(A)))
        fans.append([n,len(A),h])
    for a,b in itertools.product(range(1,101),repeat=2):
        expect=2**math.gcd(a,b)+1 if v2(a)==v2(b) else 1
        assert math.gcd(2**a+1,2**b+1)==expect
        R=lambda t:(4**t-1)//3
        assert math.gcd(R(a),R(b))==R(math.gcd(a,b))
        assert (R(b)%R(a)==0)==(b%a==0)
        assert math.gcd(4**a-1,4**b-1)==4**math.gcd(a,b)-1
    for n in range(1,101):
        assert sum(math.gcd(4**a-1,4**n-1)==3 for a in range(1,n+1))==sum(math.gcd(a,n)==1 for a in range(1,n+1))
    minima=[]
    for k in range(1,6):
        feasible=[c for c in itertools.combinations(range(1,17),k) if len({v2(a) for a in c})==k]
        best=min(sum(a+2 for a in c) for c in feasible)
        winners=[c for c in feasible if sum(a+2 for a in c)==best]
        assert best==2**k-1+2*k and winners==[tuple(2**i for i in range(k))]
        minima.append([k,best,list(winners[0])])
    families=failures=0
    for r in range(1,4):
        for weights in itertools.combinations_with_replacement(range(1,13),r):
            families+=1;fs=[factors(w) for w in weights]
            for i,f in enumerate(fs):
                private=set(f)-set().union(*(set(g) for j,g in enumerate(fs) if j!=i))
                if not private:
                    # Explicit finite-support witness, even for unit weights.
                    e=[0]*r
                    for p,power in f.items():
                        j=next(j for j,g in enumerate(fs) if j!=i and p in g)
                        e[j]=max(e[j],(power+fs[j][p]-1)//fs[j][p])
                    image=math.prod(w**power for w,power in zip(weights,e))
                    assert e[i]==0 and image%weights[i]==0
                    failures+=1
            coprime=all(math.gcd(weights[i],weights[j])==1 for i in range(r) for j in range(i))
            # Gcd/lcm preservation on different unit coordinates is already necessary.
            unit_tests=all(math.gcd(weights[i],weights[j])==1 and math.lcm(weights[i],weights[j])==weights[i]*weights[j] for i in range(r) for j in range(i))
            assert coprime==unit_tests
    assert 51*85%15==0 and math.gcd(15,51)==3
    # Normalization is not a realizability assertion: R_3=21 is a forbidden H.
    assert (4**3-1)//3==21
    def phi_value(m):
        return math.prod((p-1)*p**(e-1) for p,e in factors(m).items())
    image_totients=[]
    for n in range(1,17):
        v=phi_value(4**n-1)
        assert (v&(v-1)==0)==(n in {1,2,4,8,16})
        image_totients.append([n,v])
    assert (4**32-1)%641==0 and 640%5==0
    prime_checks=0
    for p in range(3,102,2):
        if factors(p)!={p:1}:continue
        d=next(t for t in range(1,p) if pow(4,t,p)==1)
        def valuation(m):
            e=0
            while m%p==0:e+=1;m//=p
            return e
        for n in range(1,201):
            assert ((4**n-1)%p==0)==(n%d==0)
            if n%d==0:assert valuation(4**n-1)==valuation(4**d-1)+valuation(n//d)
            prime_checks+=1
    A,B=fan(2),fan(3)
    witness=next(P for P in itertools.permutations(range(len(B)),len(A)) if all(A[i][j]==B[P[i]][P[j]] for i in range(len(A)) for j in range(len(A))))
    assert 3%2!=0 and (4**3-1)%(4**2-1)!=0
    # floor(log_4(X+1)) evaluated by exact integer thresholds, not floats.
    for X in range(4097):
        floor_log=0;power=4
        while power<=X+1:floor_log+=1;power*=4
        assert floor_log==sum(4**n-1<=X for n in range(1,8))
    out={'status':'FINITE-EXACT independent checks PASS','near_rows_a_order_H':near_rows,'fan_rows_n_order_H':fans,'arithmetic_pair_universe':'1 <= a,b <= 100','totient_n_universe':'1 <= n <= 100','finite_optima_exponents_1_to_16':minima,'scalar_weight_universe':'multisets of size 1..3 over 1..12','weight_families':families,'constructed_failed_private_coordinate_witnesses':failures,'normalization_hostile_R3':21,'producer_imports':False,'image_totients_n_1_to_16':image_totients,'independent_induced_Q2_in_Q3':list(witness),'prime_visibility_and_LTE_checks':prime_checks,'sparse_count_X_universe':'integers 0..4096'}
    Path(__file__).with_name('verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
