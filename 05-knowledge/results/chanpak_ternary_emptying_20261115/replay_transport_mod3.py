from collections import Counter
from itertools import product
base=Counter((-8,-24,17,1))
step=Counter(8+sum(i*j for i,j in zip(bs,(-16,25,25,-16,-34))) for bs in product((0,1),repeat=5))
prim=step.copy();prim.subtract(base)
def check(t):
    if not t:raise AssertionError
def mul(a,b):
    N=Counter()
    for x,i in a.items():
        for y,j in b.items():N[(x+y)%3]+=i*j
    return N
for typ,A in [('all',step),('primitive',prim)]:
    W=Counter((x%3 for x in (-8,-24,17,1)))
    for k in range(1,201):
        L=sum(W.values())
        if typ=='all':
            expect=[(4*32**(k-1)+(-1)**(k-1)*(3*(r==1)-1))//3 for r in range(3)]
        else:
            expect=[(4*28**(k-1)+3*(r==(2*k-1)%3)-1)//3 for r in range(3)]
        check([W[r] for r in range(3)]==expect and sum(expect)==L)
        W=mul(W,A)
print('exact ternary displacement coefficients: 400 width controls')
