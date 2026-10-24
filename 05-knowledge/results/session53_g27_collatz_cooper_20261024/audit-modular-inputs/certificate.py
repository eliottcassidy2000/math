"""Independent exact q-series certificate; standard library, no disabled asserts."""
import json
from math import gcd
N=96
unit=[1]+[0]*N

def require(ok, message):
    if not ok: raise ValueError(message)
def plus(a,b): return [x+y for x,y in zip(a,b)]
def scale(a,k): return [k*x for x in a]
def times(a,b):
    out=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j in range(N+1-i): out[i+j]+=x*b[j]
    return out

def power(a,n):
    out=unit[:]
    while n:
        if n&1: out=times(out,a)
        n//=2
        if n: a=times(a,a)
    return out

def inverse(a):
    require(a[0]==1,'inverse requires unit constant 1')
    b=unit[:]
    for n in range(1,N+1): b[n]=-sum(a[k]*b[n-k] for k in range(1,n+1))
    return b

def dilate(a,n): return [a[k//n] if k%n==0 else 0 for k in range(N+1)]
def chi(n): return 0 if n%7==0 else (1 if n%7 in (1,2,4) else -1)
# Theta independently from the binary quadratic form, not Eisenstein sums.
theta=[0]*(N+1)
for x in range(-2*N,2*N+1):
    for y in range(-N,N+1):
        n=x*x+x*y+2*y*y
        if n<=N: theta[n]+=1
# Bounds are exhaustive because Q=(x+y/2)^2+7y^2/4.
euler=unit[:]
for d in range(1,N+1):
    for n in range(N,d-1,-1): euler[n]-=euler[n-d]
E=[0]+power(times(euler,dilate(euler,7)),3)[:-1]
h=[0]+power(times(dilate(euler,7),inverse(euler)),4)[:-1]
t=times(h,inverse(plus(plus(unit,scale(h,13)),scale(times(h,h),49))))
require(times(t,power(theta,3))==E,'S4 exact coefficient check')
# Pole-cleared form, computed WITHOUT rational t or its denominator.
M=plus(times(power(theta,9),dilate(E,3)),scale(times(power(E,3),power(dilate(theta,3),3)),-1))
require(all(x%9==0 for x in M),'weight12 level21 form fails mod9')
raw=plus(dilate(t,3),scale(power(t,3),-1))
require(all(x%9==0 for x in raw),'direct exceptional lift check fails')
require(M==times(times(power(theta,9),power(dilate(theta,3),3)),raw),'clearing identity mismatch')
A1=[0]*(N+1); A2=[0]*(N+1)
for d in range(1,N+1):
    for n in range(d,N+1,d):
        A1[n]+=chi(d)*d*d
        A2[n]+=chi(n//d)*d*d
left=plus(plus(scale(unit,8),scale(A1,-7)),scale(A2,49))
right=plus(scale(power(theta,3),8),scale(E,-6))
require(left==right,'weight3 companion identity fails')
# Generalized Bernoulli constant for Eisenstein E_{3,1,chi}.
from fractions import Fraction as F
B3=lambda x:x**3-F(3,2)*x*x+F(1,2)*x
B3chi=49*sum((chi(a)*B3(F(a,7)) for a in range(1,8)),F(0))
require(-B3chi/6==F(-8,7),'Eisenstein constant mismatch')
# Eta cusp orders at the four cusps of X0(21), denominator c.
def eta_order(c,exponents):
    return sum((F(21*gcd(c,d)**2*r,24*gcd(c*c,21)*d) for d,r in exponents.items()),F(0))
cusps={str(c):{'E':str(eta_order(c,{1:3,7:3})),'E_at_3tau':str(eta_order(c,{3:3,21:3}))} for c in (1,3,7,21)}
require(all(F(v)>0 for item in cusps.values() for v in item.values()),'eta cusp nonholomorphy')
result={'status':'PASS','range_q_inclusive':N,'exceptional':{'weight':12,'level':21,'character':'trivial','index':32,'sturm_bound':32,'modulus':9,'coefficients_q0_to_q32':M[:33],'coefficients_divided_by_9_q0_to_q32':[x//9 for x in M[:33]],'all_residues_q0_to_q96':[x%9 for x in M]},'companion':{'weight':3,'level':7,'character':'chi_minus7','index':8,'sturm_bound':2,'eight_times_coefficients_q0_to_q2':left[:3],'E_3_1_chi_constant':str(-B3chi/6)},'eta_cusp_orders_level21':cusps,'not_a_claim':'Finite computations alone do not prove the arbitrary-prime symbolic lemmas; see REPORT.md.'}
print(json.dumps(result,sort_keys=True,indent=2))
