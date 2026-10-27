"""Independent exact reductions for THM-4812, frozen source 87d61. No assert-dependent checks."""
from fractions import Fraction
from math import factorial, log
from itertools import product
import json

def require(ok, message):
    if not ok: raise RuntimeError(message)
def f(n):
    q,r=divmod(n,3)
    return (3*5**(q-1),5**q,9*5**(q-1))[r]
def B(n,tc):
    terms=[]
    for r,p in product(range(3),range(2)):
        if 3**r*9**p >tc or 5*r+7*p>n: continue
        terms.append(3**((n-5*r-7*p)//3)*5**r*21**p)
    return max(terms,default=0)
window={v:[n for n in range(3,200) if f(n)<=v*B(n,v)] for v in (7,21)}
require(window[7]==list(range(3,13)), 'tc7 finite window')
require(window[21]==list(range(3,21))+[22], 'tc21 finite window')
# For n>=10 all admissible exceptional orbit selections fit and each term scales by 3 at n+3.
for v in (7,21):
    for n in range(10,200):
        require(B(n+3,v)==3*B(n,v), 'B recurrence')
        require(f(n+3)==5*f(n), 'f recurrence')
    cut=13 if v==7 else 23
    for n in range(cut,cut+3): require(f(n)>v*B(n,v),'tail bases')
configs=[]
for n in list(range(13,21))+[22]:
    for r,p in ((0,0),(1,0),(2,0),(0,1)):
        for c in range(n//3+1):
            a=n-3*c-5*r-7*p
            if a<0:continue
            P=3**c*5**r*21**p
            if 21*P>=f(n):configs.append(dict(n=n,k1=a,c3=c,r5=r,p7=p,P=P,index_cap=21*P//f(n)))
require(len(configs)==24,'24 configs')
require(max(c['index_cap'] for c in configs)==6,'index cap')
require(not any(c['P']%5==0 and c['index_cap']>=5 for c in configs),'index 5 possible')
# Greedy independent-set bound ends BEFORE a nonpositive factor, never resumes after it.
greedy=[]
for m in range(19,402,2):
    C=m*(m*m-1)//24
    E=3*((m*m-1)//8-1)+1
    term=Fraction(1); total=term; k=0
    while C-k*E>0:
        k+=1
        term=term*2*(C-(k-1)*E)/k
        total+=term
    require(total>23*3**((m-1)//2), f'greedy bound fails at {m}')
    greedy.append(dict(m=m,terms=k+1,bound_num=total.numerator,bound_den=total.denominator))
# For k=floor(m/20)>=20: 2^k i_k >= ((11/12)m^2)^k >= ((1100/3) k^2)^k.
# m<=20k+19 gives 23*3^((m-1)/2)<=23*3^(10k+9).
# ((1100/3) k^2 /3^10)^k increases with k>=20; test its base exactly.
require((1100*20**2)**20 > 23*3**9*(3**11)**20, 'analytic-tail base')
# Literal monotonicity sentence in source is false: fixed-k intervals have a descending margin.
def margin(m):
    k=m//20
    return k*log(.0458*m**3)-log(factorial(k))-log(23)-(m-1)/2*log(3)
require(403**60 < 3*401**60, 'expected exact source monotonicity defect')
require(margin(403)<margin(401), 'expected source monotonicity defect')
result=dict(status='PASS with proof repairs',window=window,configs=configs,greedy_min_margin=min((float(Fraction(x['bound_num'],x['bound_den']))/(23*3**((x['m']-1)//2)),x['m']) for x in greedy),greedy_checks=len(greedy),tail_exact_base=True,source_margin_401=margin(401),source_margin_403=margin(403),source_margin_419=margin(419),proof_notes=['B(n+3)=3B(n) once every allowed R5/P7 selection fits (n>=10); f(n+3)=5f(n), so three failing residue bases prove infinite exclusion.','The only cap>=5 config has n13,C3^2,P7 and P=189, so index5 cannot divide P.','Greedy sum truncates at first nonpositive factor; for m>=400 the factorial-free integer tail argument proves the claimed lower bound.'])
print(json.dumps(result,sort_keys=True,indent=2))
