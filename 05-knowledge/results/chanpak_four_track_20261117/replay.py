from base_replay import hams,validate,prime
from itertools import product
def profile(a):
    validate(a)
    h=hams(a);return h[-1],sum(x*h[len(h)-1-i] for i,x in enumerate(h))
def flip(a,i,j):
    r=list(a);r[i]^=1<<j;r[j]^=1<<i;return tuple(r)
A=(134,124,152,161,169,133,61,66)
B=(22,68,104,35,238,67,137,47)
C=(22,84,120,115,96,67,1)
a=[profile(A),profile(flip(A,1,4))]
b=[profile(B),profile(flip(B,0,2))]
h,p=profile(C)
if (a,b,h,p)!=([(485,6968),(491,6968)],[(359,5528),(383,5528)],77,1108):
    raise ValueError("module")
P=2*a[0][1]*b[0][1]*p
states=[h*x[0]*y[0] for y,x in product(b,a)]
vals=[68*P+5*z for z in states]
if vals!=[5804433777827,5804434607117,5804438259227,5804439143957]:
    raise ValueError(vals)
for p in vals:
    if not prime(p):raise ValueError(p)
    print(p,7*p)
D=vals[0]+vals[3]-vals[1]-vals[2]
if D!=55440 or vals[0]*vals[3]-vals[1]*vals[2]!=68*P*D:
    raise ValueError("inverse")
print("P2",P,"offset",68*P,"mixed difference",D)
