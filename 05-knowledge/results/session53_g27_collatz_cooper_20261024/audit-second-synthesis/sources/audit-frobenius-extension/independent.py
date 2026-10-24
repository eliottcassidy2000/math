import json
from math import comb, isqrt

def require(ok, why):
    if not ok: raise RuntimeError(why)
N=240
# Independent definition, not the producer's three-term recurrence.
F=[sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range((n+1)//2,n+1)) for n in range(N+1)]
c=[1]
for n in range(1,N+1):
    v,r=divmod(F[n]-sum(c[k]*c[n-k] for k in range(1,n)),2)
    require(r==0,'square root integrality');c.append(v)
require(c[:5]==[1,2,22,336,6006],'normalization')
ps=[p for p in range(2,isqrt(N)+1) if all(p%d for d in range(2,isqrt(p)+1))]
rows=[]
for p in ps:
    for n in range(N//(p*p)+1): require((c[p*p*n]-c[n])%(p*p)==0,'two-step')
    rows.append([p,N//(p*p)])
w=[1,2,4,3,3,6,0,0,0]
def digit(n):
    a=1
    while n:
        n,d=divmod(n,9);a=a*w[d]%9
    return a
require(all(c[n]%9==digit(n) for n in range(N+1)),'mod9 digit product')
counts=[]
hcycle=[2,1,-1,-2,-1,1]
for r in range(1,6):
    vals=[digit(n) for n in range(9**r)]
    units=sum(v%3!=0 for v in vals)
    v1=sum(v in (3,6) for v in vals)
    require(units==3**r and v1==r*3**r,'digit budget')
    hist=[vals.count(v) for v in range(9)]
    for j in range(6):
        numer=3**r+(-1)**j+2**r*hcycle[(r-j)%6]
        require(numer%6==0 and hist[pow(2,j,9)]==numer//6,'exact unit histogram')
    require(hist[3]==(r*3**r+r)//2 and hist[6]==(r*3**r-r)//2,'defect histogram')
    require(hist[0]==9**r-(r+1)*3**r,'zero histogram')
    counts.append([r,units,v1,units+v1,hist])
def mul(a,b,m):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%m
    return out
root=[1,5,1]
require(mul(root,root,9)==[1,1,0,1,1],'square root product polynomial')
p=mul(root,[1,6,0,6],9)
require(p==w[:6],'digit polynomial expansion')
require((c[9]-c[1])%27==18,'mod27 hostile')
require((c[3]-15)%9==6,'false simple swap hostile')
# U acts on all integral indices. Verify CT adjoint with both Laurent tails.
A={i:(i*i+3*i+7) for i in range(-23,24)}
B={i:(2*i*i-i+1) for i in range(-11,12)}
for p in (3,5,7):
    lhs=sum(a*B.get(-i//p,0) for i,a in A.items() if i%p==0)
    UA={i//p:a for i,a in A.items() if i%p==0}
    rhs=sum(a*B.get(-i,0) for i,a in UA.items())
    require(lhs==rhs,'Laurent CT adjoint')
    wrong=sum(a*B.get(-i,0) for i,a in UA.items() if i>=0)
    require(wrong!=lhs,'negative-index truncation hostile')
# Corrupt digit weight and normalization: acceptance must fail, including -O.
bad=w[:];bad[2]=5
require(c[2]%9!=bad[2],'digit corruption detected')
require(c[2]**2!=F[2],'coefficientwise-root hostile')
print(json.dumps({'status':'PASS','scope':'FINITE-EXACT supplemental tests, not proof of infinite statements','binomial_square_root_bound':N,'two_step_rows_p_nmax':rows,'mod9_digit_bound':N,'digit_counts_r_units_v1_nonzero':counts,'c9':c[9],'mod27_failure':18,'hostiles':['simple swap mod9','mod27 strengthening','negative-index deletion','digit mutation','coefficientwise square root']},indent=2,sort_keys=True))
