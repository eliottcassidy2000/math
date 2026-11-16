from decoder import inverse,append_word,CORE,ROOT
from itertools import permutations
from math import isqrt

def seed_profile(rows,root):
    n=len(rows); full=(1<<n)-1
    D=[[0]*n for _ in range(full+1)]
    H=[1]+[0]*full
    for s in range(1,full+1):
        for v in range(n):
            if not s>>v&1: continue
            sub=s^(1<<v)
            D[s][v]=sum(D[sub][j] for j in range(n) if rows[v]>>j&1) if sub else 1
        H[s]=sum(D[s])
    E=sum(D[s][root]*H[full^s] for s in range(1,full+1) if s.bit_count()%2==0)
    O=sum(D[s][root]*H[full^s] for s in range(1,full+1) if s.bit_count()%2)
    return H[-1],H[full^(1<<root)],E,O

def canon(rows):
    from itertools import product
    n=len(rows)
    cells=[tuple(i for i in range(n) if rows[i].bit_count()==k) for k in range(n)]
    return min(tuple(sum((rows[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
        for choice in product(*(permutations(c) for c in cells if c))
        for p in [sum(choice,())])

def seeds(limit=7):
    out={}
    layer={(0,)}
    for n in range(1,limit+1):
        for rows in sorted(layer):
            for r in range(n):
                h,g,e,o=seed_profile(rows,r)
                if o==e+g and e>0:
                    out.setdefault((2*e,g,h-2*e-g),(rows,r))
        if n<limit:
            layer={canon(tuple(row|((1-(mask>>i&1))<<n) for i,row in enumerate(rows))+(mask,))
                for rows in layer for mask in range(1<<n)}
    return out

def primality(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))

if __name__=='__main__':
    bank=seeds()
    print('bank size',len(bank),flush=True)
    for p in [30011,30013,1000000007,1000000009]:
        if not primality(p):raise ValueError(p)
        hits=[]
        from math import gcd
        for (A,g,j),b in bank.items():
            if (p-j)%gcd(A,g): continue
            if res:=list(inverse(p-j,A,g)):
                hits.append((A,g,j,b,res))
        print(p,hits[:5],len(hits),flush=True)

def transfer(profile,z):
    H,g,E,O=profile
    for bit in reversed(bin(z)[2:]) if z else ():
        if bit=='1': H,E,O=H+2*E,E+O-g,E+O
        else: E,O=E+O,E+O
    return H,g,E,O

WITNESSES=(
 (30011,(2,4,9,19,39,15),2,1561),
 (30013,(2,4,9,19,39,15),1,2363),
 (1000000007,(2,4,9,19,39,15),2,47827936),
 (1000000009,(2,12,25,49,35,7),1,30268243),
)

def literal_checks():
    from decoder import address_count,hamiltonian
    for p,rows,r,z in WITNESSES:
        h,g,e,o=seed_profile(rows,r)
        if not (primality(p) and o==e+g and hamiltonian(rows)==h
                and address_count(z,2*e,g)+h-2*e-g==p
                and transfer((h,g,e,o),z)[0]==p):
            raise ArithmeticError(p)
        # Independent subset DP verifies initial address prefixes.
        for k in range(6):
            b=z&((1<<k)-1)
            if hamiltonian(append_word(rows,r,b))!=transfer((h,g,e,o),b)[0]:
                raise ArithmeticError((p,k))
        print(p,len(rows)+z.bit_length(),(h,g,e,o),z,transfer((h,g,e,o),z))
    tri=(4,1,2)
    for z,p in zip((289,297,278,286), (721,723,725,727)):
        if hamiltonian(append_word(tri,0,z))!=p:
            raise ArithmeticError(z)
    print('prescribed prime programs and interleaved patch passed')
