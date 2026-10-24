import json

def require(x,s):
    if not x:raise RuntimeError(s)
N=8000
c=[1,2]
for n in range(1,N):
    a=(26*n*n+13*n+2)*c[n]+3*(3*n-1)*(3*n-2)*c[n-1]
    q,r=divmod(a,(n+1)**2);require(r==0,'integer recurrence');c.append(q)
weights=[1,2,4,3,3,6,0,0,0]
def prediction(n):
    product=1;defects=0;forbidden=False
    while n:
        n,d=divmod(n,9);product=product*weights[d]%9
        defects+=3<=d<=5;forbidden|=d>=6
    return product,not forbidden and defects<=1
for n,v in enumerate(c):
    pred,support=prediction(n)
    require(v%9==pred,'digit product');require((v%9!=0)==support,'one defect support')
# Exact polynomial multiplication, all arithmetic modulo9.
a=[1,5,1];b=[1,6,0,6];P=[0]*6
for i,x in enumerate(a):
    for j,y in enumerate(b):P[i+j]=(P[i+j]+x*y)%9
require(P==weights[:6],'polynomial multiplication')
require(all((c[9*n]-c[n])%9==0 for n in range(N//9+1)),'C9 fixed point')
require(all((c[27*n]-c[3*n])%9==0 for n in range(N//27+1)),'C27=C3')
# C3G=(1+3t)G(t3): coefficient of t^n on right.
for n in range(N//3+1):
    rhs=(c[n//3] if n%3==0 else 3*c[(n-1)//3] if n%3==1 else 0)
    require((c[3*n]-rhs)%9==0,'corrected companion')
# H=G/s. Q=1-26t-27t²; construct its exact integer root and H through 20.
q=[1,-26,-27]+[0]*18;s=[1]
for n in range(1,21):
    z=q[n]-sum(s[k]*s[n-k] for k in range(1,n));require(z%2==0,'s integral');s.append(z//2)
h=[1]
for n in range(1,21):h.append(c[n]-sum(s[k]*h[n-k] for k in range(1,n+1)))
require(h[:3]==[1,15,315],'H normalization')
for n in range(21):
    rhs=(c[n//3] if n%3==0 else 6*c[(n-1)//3] if n%3==1 else 0)
    require((h[n]-rhs)%9==0,'H companion')
require((c[3]-h[1])%9==6,'naive companion fails')
require((c[9]-c[1])%27!=0,'mod27 fixed point hostile')
corrupt=weights[:];corrupt[3]=6
require(corrupt[3]!=c[3]%9,'corruption detection')
print(json.dumps({'status':'FINITE-EXACT checks supporting separately stated all-index proof','through':N,'digit_weights':weights,'polynomial_product':P,'one_defect_criterion_checked':True,'corrected_C3_companion_checked_through':N//3,'H_companion_checked_through':20,'mod27_n1_failure':{'c9':c[9],'c1':c[1],'difference_mod27':(c[9]-c[1])%27},'naive_C3G_equals_H_n1_difference_mod9':(c[3]-h[1])%9},sort_keys=True,indent=2))
