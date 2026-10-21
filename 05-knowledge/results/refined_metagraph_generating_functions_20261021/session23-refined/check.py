"""Exact formal-alphabet checks; not a tournament census or atom-realizability claim."""
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
import hashlib, json, math
import sympy as sp

N=13
HMAX=81
KMAX=20
checks=0

def need(test, message):
    global checks
    checks += 1
    if not test:
        raise ValueError(message)

# name, order, h, local blue, converse name, aut-v3, aut-v5.
# K1/C3/R5/S4 signatures are actual. a,b,c,d are synthetic hostile symbols.
# Synthetic h=9 pure-blue d is NOT asserted realizable as a strong tournament.
letters=[('1',1,1,1,'1',0,0),('t',3,1,1,'t',1,0),
         ('r',5,3,3,'r',0,1),('s',4,5,3,'s',0,0),
         ('a',2,7,0,'b',1,0),('b',2,7,0,'a',1,0),
         ('c',3,9,1,'c',0,1),('d',4,9,9,'d',1,1)]
D={a[0]:a for a in letters}
zero=(0,0,0,0)
one=Counter({zero:1})

def add(*polys):
    out=Counter()
    for p in polys:
        for key,val in p.items():out[key]+=val
    return Counter({key:val for key,val in out.items() if val})

def mul(a,b):
    out=Counter()
    for x,c in a.items():
        for y,d in b.items():
            key=tuple(x[i]+y[i] for i in range(4))
            if key[0]<=N:out[key]+=c*d
    return out

def square_substitute(a):
    return Counter({tuple(2*i for i in key):val for key,val in a.items() if 2*key[0]<=N})

def mon(key):return Counter({key:1})

def weight(letter):
    _,n,h,b,opp,v3,v5=letter
    return (n,1,v3,v5)

# Closed binomial formula for neutral words, independent of enumeration.
A=Counter()
for n in range(N+1):
    for triangles in range(n//3+1):
        singles=n-3*triangles
        A[(n,singles+triangles,triangles,0)]=math.comb(singles+triangles,triangles)
P=defaultdict(Counter)
centers=[(zero,1,1)]
for letter in letters:
    name,n,h,b,opp,v3,v5=letter
    P[h][weight(letter)]+=1
    if name==opp:centers.append((weight(letter),h,b))

@lru_cache(None)
def factors(h):
    if h==1:return ((),)
    result=[]
    for j in sorted(P):
        if j>1 and h%j==0:
            result.extend((j,)+tail for tail in factors(h//j))
    return tuple(result)

F={}
for h in range(1,HMAX+1,2):
    summands=[]
    for f in factors(h):
        poly=A
        for j in f:poly=mul(mul(poly,P[j]),A)
        summands.append(poly)
    F[h]=add(*summands)

# Literal words and explicit reverse/opposition quotient, a separate path.
ordinary=defaultdict(Counter)
selfcon=defaultdict(Counter)
selfdeficit=defaultdict(Counter)
merged=defaultdict(Counter)
seen=set()
word_count=0
stack=[((),zero,1)]
while stack:
    word,key,h=stack.pop()
    word_count+=1
    ordinary[h][key]+=1
    opp=tuple(D[s][4] for s in word[::-1])
    orbit=min(word,opp)
    if word==opp:
        x=math.prod(D[s][2] for s in word[:len(word)//2])
        beta=1 if len(word)%2==0 else D[word[len(word)//2]][3]
        b=x*beta
        selfcon[h][key]+=1
        selfdeficit[h-b][key+(h,b)]+=1
        k=(h-b)//2
        if orbit not in seen:merged[k][key+(h,b)]+=1
    elif orbit not in seen:
        merged[h][key+(h,0)]+=1
    seen.add(orbit)
    for letter in letters:
        name,n,hh,b,oo,v3,v5=letter
        nxt=(key[0]+n,key[1]+1,key[2]+v3,key[3]+v5)
        if nxt[0]<=N and h*hh<=HMAX:stack.append((word+(name,),nxt,h*hh))

J={}
for h in range(1,HMAX+1,2):
    need(F[h]==ordinary[h],f'F marked h={h}')
    terms=[]
    for x in range(1,math.isqrt(h)+1,2):
        if h%(x*x):continue
        for wk,a,beta in centers:
            if x*x*a==h:terms.append(mul(square_substitute(F[x]),mon(wk)))
    J[h]=add(*terms)
    need(J[h]==selfcon[h],f'J marked h={h}')

for k in range(1,KMAX+1):
    d=2*k
    s=Counter()
    for x in range(1,d+1,2):
        if d%x or x*(x-1)>d:continue
        for wk,a,beta in centers:
            if x*(x*a-beta)==d:
                need(x*x*a<=HMAX,f'bounded h universe k={k}')
                for key,val in mul(square_substitute(F[x]),mon(wk)).items():
                    s[key+(x*x*a,x*beta)]+=val
    need(s==selfdeficit[d],f'S marked d={d}')
    result=Counter(s)
    if k%2:
        for key in F[k].keys()|J[k].keys():
            diff=F[k][key]-J[k][key]
            need(diff%2==0,'coefficient parity under converse')
            if diff:result[key+(k,0)]+=diff//2
    need(result==merged[k],f'B marked k={k}')

# An actual h=1 palindrome with a C3 pair: order6,SCC2,nu3(Aut)=2.
need(J[1][(6,2,2,0)]==1,'paired markers must double')
wrong=Counter({(2*n,c,v3,v5):v for (n,c,v3,v5),v in F[1].items() if 2*n<=N})
need(wrong[(6,2,2,0)]==0,'z-only substitution hostile rejected')
need(merged[3][(10,2,0,2,9,3)]==1,'R5 opposite-pair hostile retained')
need(merged[0][(10,2,0,2,9,9)]==0,'nonrigid pure-blue half is not harmless')

# Exact irreducibility, coprimality, and modulo-two noncancellation controls.
z=sp.symbols('z')
p=1-z-z**3
q=1-z*z-z**6
need(sp.Poly(p,z).is_irreducible,'p irreducible')
need(sp.Poly(q,z).is_irreducible,'q irreducible')
need(sp.gcd(p,q)==1,'p q coprime')
need(sp.Poly(q-p*p,z,modulus=2).is_zero,'q reduces to p squared')
for r in range(41):
    m=r//2
    L=z**(10*m)*(1+z+z**3) if r%2==0 else z**(10*m+5)
    need(not sp.Poly(L,z,modulus=2).rem(sp.Poly(q,z,modulus=2)).is_zero,
         f'mod2 noncancellation r={r}')
    if r!=1:
        for t in range(r+1):
            if 3**t*(3**t-1)<=2*3**r:need(t<=m,f'feasible half bound r={r}, t={t}')

# Synthetic rational F/J/S built independently from explicit SymPy products.
# Sparse symbols respect the required actual h=1 and h=3 alphabets, while
# higher atoms are formal. They stress composite factorizations and centers.
Puni={j:sum(val*z**key[0] for key,val in poly.items()) for j,poly in P.items()}
@lru_cache(None)
def Frat(h):
    if h==1:return 1/p
    return sp.cancel(sum(Frat(h//j)*Puni[j]/p for j in Puni if j>1 and h%j==0))

def Jrat(h):
    return sp.cancel(sum(Frat(x).subs(z,z*z)*z**wk[0]
        for x in range(1,math.isqrt(h)+1,2) if h%(x*x)==0
        for wk,a,beta in centers if x*x*a==h))

def Srat(d):
    return sp.cancel(sum(Frat(x).subs(z,z*z)*z**wk[0]
        for x in range(1,d+1,2) if d%x==0 and x*(x-1)<=d
        for wk,a,beta in centers if x*(x*a-beta)==d))

orders=[]
for r in range(8):
    h=3**r
    e=2 if r==1 else r//2+1
    # In this alphabet every attainable 3-power factors using h=3 or h=9.
    b=sp.cancel(Srat(2*h)+(Frat(h)-Jrat(h))/2)
    den=sp.fraction(b)[1]
    den=sp.Poly(den,z).monic().as_expr()
    expected=sp.Poly(p**(r+1)*q**e,z).monic().as_expr()
    need(sp.expand(den-expected)==0,f'exact degree3power denom r={r}')
    orders.append(sp.degree(den,z))

# Formula (10) coefficient check uses direct neutral words, not fitting.
for c in range(N//3+1):
    for n in range(3*c,N+1):
        need(ordinary[1][(n,n-2*c,c,0)]==math.comb(n-2*c,c),
             f'raw automorphism coefficient c={c},n={n}')

receipt={
 'status':'PASS','scope':'formal finite involutive alphabet; not a tournament census',
 'max_order':N,'h_cutoff':HMAX,'degree_cutoff':KMAX,'literal_words':word_count,
 'checks':checks,'variables':['order','SCC count','nu3 Aut','nu5 Aut','h','blue'],
 'exact_denominator_3power_r_through':7,'recurrence_orders':list(map(int,orders)),
 'mod2_top_numerator_r_through':40,
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'proof_sha256':hashlib.sha256(Path(__file__).with_name('PROOF.md').read_bytes()).hexdigest(),
 'synthetic_letters':letters
}
print(json.dumps(receipt,indent=2,sort_keys=True))
