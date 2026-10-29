"""Independent even-subset symplectic enumeration and fixed-code audit."""
from itertools import combinations, product
import json

def need(b,m):
    if not b:raise ArithmeticError(m)

def canon(x):return min(x,x^255)
def orth(u,v):return (u&v).bit_count()%2==0
def rot(x):return ((x&127)<<1&127)|((x>>6)&1)|(x&128)

def main():
    # A degree-seven polynomial is reducible only if it has a factor of
    # degree at most three. Exhaust all monic factors over F5 independently.
    f=[1,1,0,0,0,0,0,1];trials=0
    for d in [1,2,3]:
        for coeff in product(range(5),repeat=d):
            g=list(coeff)+[1];r=f[:]
            for k in range(7,d-1,-1):
                c=r[k]
                for j in range(d+1):r[k-d+j]=(r[k-d+j]-c*g[j])%5
            need(any(r[:d]),'irreducible branch polynomial');trials+=1
    need(trials==155,'all low-degree monic factors')
    V=[x for x in range(128) if x.bit_count()%2==0]
    need(len(V)==64,'quotient dimension6')
    spaces=set()
    for a,b,c in combinations(V[1:],3):
        if c==a^b or not all(orth(u,v) for u,v in [(a,b),(a,c),(b,c)]):continue
        L=frozenset([0,a,b,c,a^b,a^c,b^c,a^b^c])
        need(len(L)==8,'dimension3');spaces.add(L)
    need(len(spaces)==135,'all lagrangians')
    nonmatching=[]
    for L in spaces:
        C=set(L)|{x^255 for x in L}
        if any(x.bit_count()==2 for x in C):continue
        need(all(x.bit_count()%4==0 for x in C),'doubly even')
        triples=[x&127 for x in C if x&128 and x.bit_count()==4]
        need(len(triples)==7,'seven Fano lines')
        for i,j in combinations(range(7),2):
            need(sum(bool(t>>i&1 and t>>j&1) for t in triples)==1,'Fano pair incidence')
        nonmatching.append(L)
    fixed=[L for L in spaces if frozenset(canon(rot(x)) for x in L)==L]
    need(len(nonmatching)==30 and len(fixed)==2,'30 nonmatching and2fixed')
    need(all(L in nonmatching for L in fixed),'fixed kernels nonmatching')
    # Affine permutations x->ax+b in F7 preserving a fixed Fano line set.
    L=fixed[0];C=set(L)|{x^255 for x in L}
    lines={x&127 for x in C if x&128 and x.bit_count()==4}
    stabilizers=[]
    for a in range(1,7):
        for b in range(7):
            moved={sum(1<<((a*i+b)%7) for i in range(7) if t>>i&1) for t in lines}
            if moved==lines:stabilizers.append((a,b))
    need(len(stabilizers)==21 and {a for a,b in stabilizers}=={1,2,4},'normalizer21')
    print(json.dumps({'lagrangians':135,'matching':105,'extended_hamming':30,
        'cycle_1_7_fixed':2,'normalizer_size':21,'multipliers':[1,2,4],'irreducibility_factor_trials':trials},sort_keys=True))
if __name__=='__main__':main()
