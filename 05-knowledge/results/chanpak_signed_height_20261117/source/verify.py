from check import profile,blow,a,b
from math import comb,factorial
def ck(x):
    if not x:raise ArithmeticError('control failed')
def add(a,q,t):
    n=len(a);A=((1<<t)-1)<<n
    return tuple(r|(0 if q>>i&1 else A) for i,r in enumerate(a))+tuple(q|(A&~((1<<(n+i+1))-1)) for i in range(t))
def delta(x,y):
    return tuple(w-v for v,w in zip(profile(x),profile(y)))
def C(m,r):return sum((-1)**i*comb(r,i)*(r-i)**m for i in range(r+1))
for j,expected in {0:(10,-207,-170),1:(10,4,0),2:(10,4,0,0),3:(10,6,0,0),4:(10,118,54),7:(10,6,0)}.items():
    old=[]
    for r in range(1,min(a[j].bit_count(),7-a[j].bit_count())+2):
        n=delta(blow(a,j,r),blow(b,j,r))[0]-sum(d*C(r,i) for i,d in enumerate(old,1))
        ck(n%factorial(r)==0);old.append(n//factorial(r))
    ck(tuple(old)==expected)
for i,j,W,E,base in [(12,31,(7999,157792),(210,7104),1807703),(43,75,(8021,157952),(-2710,-40216),1809661)]:
    x,y=[add(add(z,i,1),j,1) for z in (a,b)]
    ck(profile(x)==W and delta(x,y)==(10,-8))
    ck(delta(*[add(add(z,i,2),j,2) for z in (a,b)])==E)
    ck(9*W[0]+11*W[1]==base)
z=(102,260,88,851,131,286,50,367,597,247)
w=list(z);w[4]^=1<<5;w[5]^=1<<4;w=tuple(w)
ck(profile(z)==(5781,119040) and profile(w)==(5611,119056))
ck(profile(add(z,0,3))==(5781,952320) and profile(add(w,0,3))==(5611,952448))
def F(t):
    if not t:return 1
    u=t//2;return F(u)+(0 if t%2 else u+u%2)
u=4160;k=18;h=1;R=64
rev=int(f'{u:018b}'[::-1],2)
ck(F(u)==4129 and rev==2080)
for _ in range(9):
    ck(R==63*h+1)
    a1,b1=3*h+3*R,3*h+4*R
    ck((a1,b1)==(192*h+3,255*h+4) and -170*a1+128*b1==2)
    ck(a1*5781+b1*952320==243951552*h+3826623)
    h,R=F(u)*h+(u+1)*R,(2**k-rev)*h+2**k*R
    ck((h-4161)%266272==0)
print('PASS source profiles, uniform-partition certificates, and independent m6 reader')
