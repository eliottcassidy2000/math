# The Fermat line under Collatz; Pascal mod 2 = odd constructible; Pepin as an affine-permutation cycle type; Paley graphs of Fermat primes.
import json, sys
from math import comb, gcd
from itertools import combinations
from sympy import n_order, isprime, factorint, totient
def v2(x): return (x & -x).bit_length() - 1
def syracuse_orbit(x, maxsteps=10**6):
    """Odd x -> list of (value, valuation) along the Syracuse map until reaching 1."""
    out=[]; steps=0
    while x != 1 and steps < maxsteps:
        y=3*x+1; v=v2(y); out.append((x,v)); x=y>>v; steps+=1
    return out
# (1) Fermat line: F_n = 2^(2^n)+1, n>=2: valuations 2 repeated 2^(n-1)-1 times, then 3, landing on (3^(2^(n-1))+1)/2.
res={'fermat_line':{}}
for n in range(1,6):
    F=2**(2**n)+1
    orb=syracuse_orbit(F, maxsteps=200000)
    vals=[v for _,v in orb]
    L=2**(n-1)-1
    pred_ok = (n==1) or (vals[:L]==[2]*L and vals[L]==3 and orb[L+1][0]==(3**(2**(n-1))+1)//2) if n>=2 else None
    # stopping time (first value below F) and total stopping time (odd steps)
    below=next((i for i,(x,_) in enumerate(orb) if x<F), None)
    res['fermat_line'][f'F_{n}={F if F<10**7 else str(F)[:12]+"..."}']={'odd_steps_to_1':len(orb) if orb and orb[-1][0]!=F else None,'run_prediction_ok':pred_ok,
        'first_valuations':vals[:2**(n-1)+3],'first_odd_step_below_start':below,'reached_1':orb[-1][0]!=1 and (3*orb[-1][0]+1)>>v2(3*orb[-1][0]+1)==1 if orb else True}
# general 2^m+1 run law: x = 1 mod 2^(2L+1) iff L consecutive valuation-2 steps
def run_len(x):
    k=0
    while k<64:
        y=3*x+1
        if v2(y)!=2: return k
        x=y>>2; k+=1
    return k
res['run_law_check']=all(run_len(x) == max(L for L in range(0,20) if (x-1)%(2**(2*L+1))==0) for x in range(3,20000,2))
res['two_power_plus_one']={m: {'run_v2':run_len(2**m+1),'predicted':(m-1)//2,'odd_steps_to_1':len(syracuse_orbit(2**m+1,10**6))} for m in range(2,41)}
# (2) odd constructible numbers (products of distinct Fermat primes) = rows of Pascal mod 2 read in binary
fermat=[3,5,17,257,65537]
prods=sorted({ (lambda S: eval('*'.join(map(str,S)) if S else '1'))(S) for r in range(0,6) for S in combinations(fermat,r)})
rows=sorted(sum(1<<k for k in range(s+1) if comb(s,k)%2==1) for s in range(32))
res['pascal_mod2_equals_fermat_products']= prods==rows
res['product_all_five']=3*5*17*257*65537
res['constructible_up_to_40']=[n for n in range(3,41) if all(p in fermat for p in factorint(n) if p!=2) and all(e==1 for p,e in factorint(n).items() if p!=2)]
# (3) Pepin as cycle type of x -> 3x+1 on Z/F_n: fixed point -1/2, one big cycle iff 3 is a primitive root iff F_n prime (n>=1)
res['affine_cycle_type']={}
for n in range(1,6):
    F=2**(2**n)+1
    o=n_order(3,F)
    res['affine_cycle_type'][f'n={n}']={'prime':isprime(F),'order_of_3':o,'phi':int(totient(F)),'single_cycle_on_units':o==F-1,
        'pepin_residue_is_minus_one': pow(3,(F-1)//2,F)==F-1}
# (4) Paley graphs on Fermat primes 5 and 17: eigenvalues (-1±sqrt p)/2, clique/independence numbers (Ramsey)
import numpy as np
for p in (5,17):
    sq={(i*i)%p for i in range(1,p)}
    A=np.array([[1 if (j-i)%p in sq else 0 for j in range(p)] for i in range(p)])
    ev=sorted(set(np.round(np.linalg.eigvalsh(A),6)))
    # clique number by brute force
    best=0
    for r in range(1,p+1):
        if any(all(A[a][b] for a,b in combinations(S,2)) for S in combinations(range(p),r)): best=r
        else: break
    res[f'paley_{p}']={'eigenvalues':ev,'predicted':[(-1-p**0.5)/2,(-1+p**0.5)/2,(p-1)/2],'clique_number':best,'self_complementary_degree':(p-1)//2}
print(json.dumps(res, indent=0, default=str))
