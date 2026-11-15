from collections import Counter
from fractions import Fraction
from itertools import product

def need(x):
    if not x: raise AssertionError("replay")
def p(n):
    s=0
    while n:
        n,b=divmod(n,2)
        if not b:s+=n+(n&1)
    return s
def r(n,k):return int(f'{n:0{k}b}'[::-1],2)
def F(x,y):return 2855+5760*x*y+4034*x+3914*y+45*(p(x)+p(y))
A=1<<24
U,V=11418043,11767483
need(128*V+r(U,24)==(A-1)*4079//45)
need(128*U+r(V,24)==(A-1)*3959//45)
need(3914*V-4034*U+45*(p(V)-p(U))==0)
for k in range(1,30):
    y=V*((A**k-1)//(A-1))
    H=F(0,y)
    for j in range(k+1):
        x=U*((A**j-1)//(A-1))
        z=V*((A**(k-j)-1)//(A-1))
        need(F(x,z)==H)
print('two example counts:',F(0,V),F(0,V*(A+1)))
base=Counter({2:1,-16:1,34:1,-48:1})
step=Counter(16+sum(a*w for w,a in zip((-32,50,50,-32,-68),bs))
             for bs in product((0,1),repeat=5))
z=base
f=[]
for k in range(1,9):
    f.append(z[0])
    n=Counter()
    for u,i in z.items():
        for v,j in step.items():n[u+v]+=i*j
    z=n
need(f==[0,10,0,5250,0,3653650,0,2834954850])
e=sum(Fraction(v,4096**k) for k,v in enumerate(f,1))
tail=Fraction(4*32**8,4096**9)/(1-Fraction(32,4096))
lo,hi=e/(1+e),(e+tail)/(1+e+tail)
need(Fraction('0.0000005960647449810')<lo<hi<Fraction('0.0000005960647449811'))
rows=(4,1,130,7,143,159,63,75)
for v in range(8):
    seen=1<<v
    for _ in range(8):
        for i in range(8):
            if seen>>i&1:seen|=rows[i]
    need(seen==255)
print('zero counts through width 96:',f)
print('first-return mass interval:',float(lo),float(hi))
