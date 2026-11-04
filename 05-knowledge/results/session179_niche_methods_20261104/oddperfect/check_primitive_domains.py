"""Standard-library exact audit. No floating point, external packages, or assertions."""
from itertools import combinations, product
from math import isqrt
import json
from pathlib import Path

def require(x, message):
    if not x: raise RuntimeError(message)
def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def order(p,r):
    t=p%r; d=1
    while t!=1:
        t=t*p%r; d+=1
    return d

def divisors(n): return [d for d in range(1,n+1) if n%d==0]
def role(p,q,n): return n%4==2 if p==q else n>=3 and n%2==1

def domains(S,q):
    result={}
    for p in S:
        D={order(p,r) for r in S if r!=p}
        result[p]=sorted(n for n in D if n>1 and role(p,q,n) and all(d in D for d in divisors(n) if d>1))
    return result

def sigma(p,n): return (p**n-1)//(p-1)
def strip(v,S):
    row={}
    for r in [2]+list(S):
        a=0
        while v%r==0:
            v//=r; a+=1
        row[r]=a
    return v,row

def exact_rows(S,q):
    rows={}
    for p,ns in domains(S,q).items():
        rows[p]=[]
        for n in ns:
            rem,w=strip(sigma(p,n),S)
            if rem==1 and w[2]==int(p==q): rows[p].append((n-1,w))
    return rows

def decide(S,q):
    R=exact_rows(S,q); found=[]; tested=0
    for chosen in product(*(R[p] for p in S)):
        tested+=1
        if all(sum(w[r] for _,w in chosen)==chosen[j][0] for j,r in enumerate(S)):
            N=1
            sig=1
            for p,(a,w) in zip(S,chosen): N*=p**a; sig*=sigma(p,a+1)
            require(sig==2*N,'false positive exact decision')
            found.append([a for a,_ in chosen])
    return tested,found

def strong(S,q):
    E={p:{r for r in S if r!=p and (order(p,r)%4!=0 if p==q else order(p,r)%2==1)} for p in S}
    for p in S:
        seen={p}; todo=[p]
        while todo:
            for r in E[todo.pop()]-seen: seen.add(r); todo.append(r)
        if seen!=set(S): return False
    return True

def main():
    primes=[p for p in range(3,60,2) if prime(p)]
    checked=0; primitive_levels=0
    # Independent primitive-factor test by trial-factorization of manageable values.
    for p in primes[:7]:
        for qrole in [False,True]:
            if qrole and p%4!=1: continue
            for n in range(2,20):
                if not (n%4==2 if qrole else n%2==1): continue
                for d in divisors(n)[1:]:
                    # Detect a primitive factor without factoring a potentially huge cofactor:
                    # strip all earlier cyclotomic-prime support from p^d-1.
                    v=p**d-1
                    from math import gcd
                    for k in range(1,d):
                        g=gcd(v,p**k-1)
                        while g>1:
                            v//=g; g=gcd(v,p**k-1)
                    require(v>1, f'primitive factor absent {p,d}')
                    primitive_levels+=1
    # Direct sigma-closed rows, independently brute-force lengths, must survive domains.
    closed=0
    for size in [2,3,4]:
        for S in combinations(primes[:8],size):
            for q in S:
                if q%4!=1: continue
                D=domains(S,q); checked+=1
                for p in S:
                    require(len(D[p])<=size-1,'domain count')
                    for n in range(2,35):
                        if not role(p,q,n): continue
                        rem,w=strip(sigma(p,n),S)
                        if rem==1:
                            closed+=1
                            require(n in D[p],f'closed row lost {S,q,p,n}')
                            require(n<=max(S)-1,'max-prime bound')
                decide(S,q)
    S=(3,5,7,11,13)
    concrete={}
    for q in [5,13]:
        require(strong(S,q),'maximal graph must pass')
        D=domains(S,q)
        require(D[7]==[],'7 must be primitive-dead')
        concrete[str(q)]={'domains':D,'maximal_graph_strong':True,'order_labels_from_7':{r:order(7,r) for r in S if r!=7}}
    # Top-level witness does not ensure divisor completeness.
    require(order(5,19)==9,'top-level witness')
    require(order(5,31)==3,'missing lower witness')
    require(sigma(5,9)==31*19*829,'factorization')
    # Completeness does not ensure row closure.
    require(domains((3,5,7,11,13),13)[5]==[5],'complete candidate')
    require(sigma(5,5)==11*71,'outside factor')
    # Exception controls demonstrate why role hypotheses matter.
    require(sigma(7,2)==8,'Mersenne d2 exception control')
    require(sigma(2,6)==63,'base2 d6 exception control')
    # Verify a multi-divisor example, with two primitive factors at one level.
    require(order(5,829)==9 and prime(829),'second primitive divisor')
    output={'support_euler_cases':checked,'direct_closed_rows_preserved':closed,'primitive_divisor_levels_checked':primitive_levels,'concrete_support':concrete,'lower_divisor_counterexample':{'p':5,'n':9,'support':[5,13,19],'missing_order':3},'all_checks':'PASS'}
    Path(__file__).with_name('RESULTS.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
if __name__=='__main__': main()
