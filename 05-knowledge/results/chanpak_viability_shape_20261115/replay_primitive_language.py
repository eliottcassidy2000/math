from collections import Counter
from itertools import product
from fractions import Fraction
base=Counter({-16:1,-48:1,34:1,2:1})
step=Counter(16+sum(b*u for u,b in zip((-32,50,50,-32,-68),w))
             for w in product((0,1),repeat=5))
step.subtract(base)
def require(x):
    if not x:raise AssertionError("primitive GF")
require(min(step.values())>=0 and sum(step.values())==28)
require(step[-16]>0 and step[16]>0)
w=base
s=Fraction(0)
coefs=[]
for k in range(1,9):
    coefs.append(w[2])
    require(sum(w.values())==4*28**(k-1))
    s+=Fraction(w[2],4096**k)
    n=Counter()
    for a,i in w.items():
        for b,j in step.items():
            n[a+b]+=i*j
    w=n
require(coefs==[1,0,178,0,76140,0,40625112,0])
print('primitive +2 multiplicities at widths 12k:',coefs)
print('code density, 8-width truncation:',float(s))
print('tail ceiling:',float(Fraction(4*28**8,4096**9)/(1-Fraction(28,4096))))
