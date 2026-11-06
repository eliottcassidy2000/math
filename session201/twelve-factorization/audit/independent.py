"""Independent physical-digit, matrix, divisor, and bounded-fiber audit.
No author-code imports. All comparisons remain active under python -O.
"""
from pathlib import Path
from fractions import Fraction
from math import isqrt
from functools import lru_cache
import hashlib,json
ROOT=Path(__file__).resolve().parent

def require(ok,why):
    if not ok:raise RuntimeError(why)
T=(1,0,2,1)
M={
'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))}
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
def mul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4))for j in range(4))for i in range(4))
def act(c,v):return tuple(sum(x*y for x,y in zip(row,v))for row in M[c])
def word(w,v=T):
    for c in w:v=act(c,v)
    return v
def mat(w):
    a=tuple(tuple(int(i==j)for j in range(4))for i in range(4))
    for c in w:a=mul(M[c],a)
    return mul(Q,a)
def binary_word(u):return ''.join('O' if x=='1' else 'Z'for x in bin(u-1)[3:])
@lru_cache(None)
def g(u):
    require(u>=1,'g index')
    return 1 if u==1 else word(binary_word(u))[0]
def tau(u):return(g(u),u%2,u,1)
def value(a,u,v):return sum(x*a[i][j]*y for i,x in enumerate(tau(u))for j,y in enumerate(tau(v)))
menus=json.loads((ROOT/'dependencies/MENUS.json').read_text())['5']
forms={r['word']:mat(r['word'])for r in menus}
require(all(forms[r['word']]==tuple(map(tuple,r['matrix']))for r in menus),'all twelve physical word matrices')
# Exact formal homogeneous-matrix classification, including every 2x2 minor.
expected_dets={'J':1,'JJ':0,'JJJ':-2,'JOZJ':3,'JJZJ':0,'JJJJ':3,'JOOZJ':5,'JJZZJ':0,'JOZJJ':-7,'JJZJJ':0,'JJJJJ':-1}
dets={};products=[];corrections={}
for w,a in forms.items():
    if not w:continue
    require(all(a[3][i]==a[i][3]==0 for i in range(4)),'nonempty homogeneous coordinates')
    z=[list(row[:3])for row in a[:3]]
    require(z[0][0]!=0,'rank-one pivot')
    h=Fraction(z[1][0]*z[0][1],z[0][0]);correction=z[1][1]-h;z[1][1]=h
    minors=[z[i][j]*z[k][l]-z[i][l]*z[k][j]for i in range(3)for k in range(i+1,3)for j in range(3)for l in range(j+1,3)]
    if all(x==0 for x in minors):products.append(w);corrections[w]=str(correction)
    aa,cc,dd,bb=a[0][0],a[0][2],a[2][0],a[2][2]
    dets[w]=aa*bb-cc*dd
    require([[aa*aa,aa*cc],[dd*aa,dd*cc+dets[w]]]==[[aa*aa,aa*cc],[aa*dd,aa*bb]],'norm-style coefficient identity')
require(dets==expected_dets,'all stated core determinants')
require(set(products)=={'JJZJ','JJZZJ'},'only two homogeneous products after parity correction')
for c,w in((3,'JJZJ'),(5,'JJZZJ')):
    vec=(c,-1,c,0)
    require(all(c*forms[w][i][j]==vec[i]*vec[j]+((c-1)if i==j==1 else 0)for i in range(4)for j in range(4)),'formal factorization identity')
# Derive g through literal binary source operations, not the author's recurrence.
MAX=65536;fibers={}
for u in range(1,MAX+1):
    a=g(u)
    require(a%2==1 and 1<=a<=u,'positive odd g')
    require((a==1)==(u&(u-1)==0),'powers of two unit fiber')
    if u>=2:
        require(word(binary_word(u))==tau(u),'full actual unit source profile')
        fibers.setdefault(u+a,[]).append(u)
    if u%2==0:require(a==g(u//2),'even recurrence')
    elif u>1:
        h=(u+1)//2
        require(a==g(h)+h-h%2,'odd recurrence')
        require(2*a>=u+3,'odd-core lower bound')
        require((2*a==u+3)==((u+1)&u==0),'lower equality family')
        require((a==u)==((u-1)&(u-2)==0),'upper equality family')

def bounded_fiber(x):
    if x%2:
        candidates=range(max(2,(2*x+2)//3),x)
        return tuple(u for u in candidates if u%2==0 and u+g(u)==x)
    candidates=range(max(3,(x+1)//2),(2*x-3)//3+1)
    return tuple(u for u in candidates if u%2 and u+g(u)==x)
for x in range(1,4097):require(bounded_fiber(x)==tuple(fibers.get(x,[])),'source-sum inverse against unrestricted digit bank')

def divisors(t):
    for d in range(1,isqrt(t)+1):
        if t%d==0:
            yield d,t//d
            if d*d!=t:yield t//d,d

def factor_inverse(n,c):
    # A single scaled-factor algorithm across four parity sheets, not author branches.
    out=set()
    for p in(0,1):
        for q in(0,1):
            for d,e in divisors(c*n-(c-1)*p*q):
                if(d+p)%c or(e+q)%c:continue
                x,y=(d+p)//c,(e+q)//c
                if x%2!=1-p or y%2!=1-q:continue
                require(x<=MAX and y<=MAX,'finite inverse control dictionary complete')
                out.update((u,v)for u in fibers.get(x,[])for v in fibers.get(y,[])if u%2==p and v%2==q)
    return out
CAP=5000;inverse_counts={}
for c,w in((3,'JJZJ'),(5,'JJZZJ')):
    direct={}
    # c>=3, x,y>=3 imply cxy-qx-py+pq >= uv. Hence this domain is complete.
    for u in range(2,CAP//2+1):
        for v in range(2,CAP//u+1):
            n=value(forms[w],u,v)
            require(n>=u*v,'strong c-form product bound control')
            if n<=CAP:direct.setdefault(n,set()).add((u,v))
    for n in range(1,CAP+1):require(factor_inverse(n,c)==direct.get(n,set()),'complete uniform divisor inverse at '+str((c,n)))
    inverse_counts[w]=dict(targets=CAP,represented_targets=len(direct),ordered_witnesses=sum(map(len,direct.values())))
# Sieve all potential prime divisors of the hostile triple.
ns=[2135743,6407227,10678711];limit=isqrt(max(ns));sieve=bytearray(b'\1')*(limit+1);sieve[0:2]=b'\0\0'
for d in range(2,isqrt(limit)+1):
    if sieve[d]:sieve[d*d::d]=b'\0'*(((limit-d*d)//d)+1)
primes=[p for p in range(2,limit+1)if sieve[p]]
for n in ns:require(n>=2 and all(n%p for p in primes if p*p<=n),'prime hostile sieve')
require(ns[1]==3*ns[0]-2 and ns[2]==5*ns[0]-4,'exact hostile triple')
require(value(forms['J'],7,191317)==ns[0],'positive J witness')
hostile_word=binary_word(191317)+'J'+binary_word(7)[::-1].translate(str.maketrans('OZ','ZO'))+'J'
s=word(hostile_word[:-1]);require(s[0]-s[1]+s[2]==ns[0],'physical hostile full profile')
# Actual dyadic realization for every fixed-core case, including both singleton cores.
height_checks=0
for a in range(1,32,2):
    for b in range(1,32,2):
        for k in range(13):
            actual={value(forms['J'],a<<r,b<<(k-r))for r in range(k+1)if(a<<r)>=2 and(b<<(k-r))>=2}
            admissible=(a>=3 and b>=3)if k==0 else k>=int(a==1)+int(b==1)
            wanted=g(a)*g(b)+(a*b-1 if k==0 else (1<<k)*a*b)
            require(actual==({wanted}if admissible else set()),'dyadic collapse and converse boundary')
            height_checks+=1
# Verify the inherited global uv/2 lower bound exactly on positive-coordinate rays.
# Columns are rays (1,0,1,0),(0,0,1,0),(2,1,2,1),(1,0,2,1).
P=((1,0,2,1),(0,0,1,0),(1,1,2,2),(0,0,1,1));U=P[2]
Pt=tuple(zip(*P));bounds={}
for w,f in forms.items():
    C=mul(mul(Pt,f),P);b=min(Fraction(4*C[i][j],U[i]*U[j])for i in range(4)for j in range(4))
    require(b>=2,'physical product lower coefficient');bounds[w]=str(b)
# Nonnegative positive-coordinate source rays plus exact formula give all-order bound.
for u in range(2,4097):
    D,E,Uu,d=tau(u);require(min(D-E-d,Uu-D+E-d,E,d-E)>=0,'unit source cone coordinates')
# Finite controls of the counting argument, without interpreting them as tail coverage.
counting=[]
for N in(16,64,256,1024):
    L=(2*N).bit_length()-1
    for G in(1,3,5,9):
        values=set();jvals=set();cores={u for u in range(1,2*G+1,2)if g(u)<=G}
        require(len(cores)<=G,'bounded number of actual odd cores')
        for a in cores:
            for b in cores:
                for r in range(L+1):
                    for ss in range(L-r+1):
                        u,v=a<<r,b<<ss
                        if min(u,v)<2 or u*v>2*N:continue
                        for w,f in forms.items():
                            n=value(f,u,v)
                            if n<=N:values.add(n)
                            if w=='J' and n<=N:jvals.add(n)
        require(len(values)<=12*G*G*(L+1)*(L+2)//2,'twelve-form counting bound')
        require(len(jvals)<=G*G*(L+1),'J collapsed counting bound')
        counting.append(dict(N=N,G=G,values=len(values),J_values=len(jvals)))
print(json.dumps(dict(status='PASS',matrix_count=len(forms),formal_core_determinants=dets,
 full_factorization_templates=products,pq_corrections=corrections,physical_g_controls=MAX,
 exact_source_sum_fibers=4096,complete_divisor_controls=inverse_counts,
 prime_triple=ns,prime_divisor_sieve_limit=limit,prime_divisor_count=len(primes),
 hostile_J_indices=[7,191317],hostile_source_profile=s,hostile_graph_word=hostile_word,
 dyadic_converse_checks=height_checks,positive_ray_product_coefficients=bounds,bounded_fiber_controls=counting,
 scope='Exact scoped lemmas and finite controls only; uniform twelve-template occupancy remains OPEN.',
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
