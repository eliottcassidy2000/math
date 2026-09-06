"""Exact modular response rank reduction; no producer imports or asserts."""
from fractions import Fraction as F
from math import factorial, prod
from pathlib import Path
import json
import sys
sys.stdout.reconfigure(newline='\n')
gates=0
def need(ok,label):
    global gates
    gates+=1
    if not ok: raise ArithmeticError(label)
def fall(n,k): return prod(n-i for i in range(k))
def row(h,z):
    return [F(factorial(h)*factorial(z+2*h),factorial(j)*factorial(z+2*j)*factorial(3*h-3*j)) for j in range(h+1)]
def endpoint_row(h,S):
    return [F(fall(h,h-j)*fall(S,2*(h-j)),factorial(3*(h-j))) for j in range(h+1)]
def reduced(q,p):
    need(q.denominator%p!=0,'integral coefficient')
    return q.numerator*pow(q.denominator,-1,p)%p
def trim(a):
    while a and not a[-1]: a.pop()
    return a
def rem(a,b,p):
    a=[v%p for v in a]
    while len(a)>=len(b):
        c=a[-1]*pow(b[-1],-1,p)%p
        d=len(a)-len(b)
        for j,v in enumerate(b): a[d+j]=(a[d+j]-c*v)%p
        trim(a)
    return a
def determinant(a,p):
    a=[r[:] for r in a]; out=1
    for j in range(len(a)):
        k=next((k for k in range(j,len(a)) if a[k][j]%p),None)
        if k is None:return 0
        if k!=j:a[k],a[j]=a[j],a[k];out=-out
        out=out*a[j][j]%p; inv=pow(a[j][j],-1,p)
        for k in range(j+1,len(a)):
            c=a[k][j]*inv%p
            for l in range(j,len(a)):a[k][l]=(a[k][l]-c*a[j][l])%p
    return out%p
def norm(b,a,p):
    n=len(b)-1
    cols=[rem([0]*j+a,b,p)+[0]*n for j in range(n)]
    return determinant([[cols[j][i] for j in range(n)] for i in range(n)],p)
def isprime(p):return p>=2 and all(p%d for d in range(2,int(p**.5)+1))
def evalrow(a,t):return sum(c*t**j for j,c in enumerate(a))

records=[]
for H in range(1,19):
    pp=[p for p in range(6*H+1,6*H+45) if isprime(p)][:2]
    need(len(pp)==2,'two primes in declared finite selection')
    for m in range(min(4,2*H-1)+1):
        for p in pp:
            z=(p+2*m+1)//2-2*H
            need(z>=1 and p>6*H and (2*z+4*H-2*m-1)%p==0,'rank-drop hypotheses')
            phi,psi=row(H,z),row(2*H,2*z)
            A=[reduced(q,p) for q in phi]; B=[reduced(q,p) for q in psi]
            R=[reduced(F(fall(2*H,m-j)*fall(2*m+1,2*(m-j)),factorial(3*(m-j))),p) for j in range(m+1)]
            need(all(A),'every first-row coefficient unit')
            need(B==[0]*(2*H-m)+R,'entire second-row factorization')
            need(A==[reduced(q,p) for q in endpoint_row(H,F(2*m+1,2))],'half-integral endpoint reduction')
            actual=(-1)**H*norm(A,B,p)%p
            small=norm(R,A,p)
            predicted=(-1)**H*pow(A[0],2*H-m,p)*small%p
            need(actual==predicted,'full determinant equals rank-m consumer')
            records.append(dict(H=H,m=m,z=z,p=p,norm=actual,residual_norm=small))
            # Periodicity in z is an actual rational coefficient check.
            if m==1:
                need([reduced(q,p) for q in row(H,z+p)]==A,'first-row periodic lift')
                need([reduced(q,p) for q in row(2*H,2*(z+p))]==B,'second-row periodic lift')
                need(R==[2*H%p,1],'rank-one residual')
                need(small==evalrow(A,-2*H)%p,'rank-one evaluation consumer')

linear=[]
for H in range(1,81):
    terms=[F(fall(H,k)*fall(F(3,2),2*k),factorial(3*k)*(2*H)**k) for k in range(H+1)]
    need(terms[0]==1 and terms[1]==F(1,16),'first two alternating terms')
    for k in range(1,H):
        need(0<terms[k+1]<terms[k]/2,'uniform decreasing positive tail')
    value=sum((-1)**k*a for k,a in enumerate(terms))
    need(F(15,16)<=value<1,'nonzero normalized half-endpoint evaluation')
    A=evalrow(endpoint_row(H,F(3,2)),-2*H)
    need(A==(-2*H)**H*value,'literal rational evaluation')
    linear.append(dict(H=H,numerator=A.numerator,denominator=A.denominator))

# A dropped residual-unit guard fails even at height2.
H,p,z=2,9601,4798
need(isprime(p),'exceptional prime is prime')
A=[reduced(q,p) for q in row(H,z)]
B=[reduced(q,p) for q in row(2*H,2*z)]
need((2*z+4*H-3)%p==0 and p>6*H,'hostile obeys rank-one endpoint hypotheses')
need(evalrow(endpoint_row(2,F(3,2)),-4)==F(9601,640),'exact exceptional rational numerator')
need(B==[0,0,0,4,1] and evalrow(A,-4)%p==0,'both rows share t+4 in hostile residue')
need(norm(A,B,p)==0,'rank-drop is not automatic modular nonvanishing')
# An independent prime separates these same rational rows.
p2=28807
need(isprime(p2) and p2>6*H,'separate positive-prime control')
C=[reduced(q,p2) for q in row(H,z)];D=[reduced(q,p2) for q in row(2*H,2*z)]
hostile_norm=(-1)**H*norm(C,D,p2)%p2
need(hostile_norm!=0,'bad endpoint reduction does not imply rational cancellation')

# An explicit new high-height carried boundary absent from the old endpoint test.
new=next(r for r in records if r['H']>=7 and r['m']==1 and r['z']%2==1 and r['norm'])
h=new['H']+(new['z']-1)//2;r=(new['z']-1)//2
need(r>=1 and (4*h-1)%new['p']==0,'actual carried-boundary dictionary')
old=[]
M=4*h+1
for p in range(2,M+1):
    if not isprime(p) or M%p or p<=4*new['H']:continue
    n=M;e=0
    while n%p==0:n//=p;e+=1
    if e>6*new['H']//p:old.append(p)
need(not old,'chosen example is not supplied by old prime-divisor condition')

certificate=dict(status='FINITE-EXACT controls for proved rank-drop identity',universe='H=1..18; m=0..min(4,2H-1); first two primes in (6H,6H+44]; half-endpoint dominance H=1..80',records=records,linear_evaluations=linear,hostile=dict(H=2,z=4798,p=9601,other_prime=p2,other_norm=hostile_norm),new_carried=dict(h=h,r=r,**new))
path=Path(__file__).with_name('continuing7_20260906_norm_rank_drop_certificate.json')
path.write_text(json.dumps(certificate,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
print('Rank-drop cases:',len(records))
print('Half-endpoint evaluation controls:',len(linear))
print('Hostile rank-one endpoint:',certificate['hostile'])
print('New carried boundary:',certificate['new_carried'])
print('Always-active gates:',gates)
