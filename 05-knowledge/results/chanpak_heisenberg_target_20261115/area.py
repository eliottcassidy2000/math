from itertools import product
from math import comb
def invariant(w):
    a=w.count('A')
    b=len(w)-a
    c=sum(x=='A' and y=='B' for i,x in enumerate(w) for y in w[i+1:])
    return a,b,c
def cycleprofile(w):
    n=len(w)+1
    h=1
    q=2**n
    for i,x in enumerate(w):
        for j in range(i+1,len(w)):
            if x!='A' or w[j]!='B':continue
            g=j-i-1
            for k in range(0,g+1,2):
                h+=2*comb(g,k);q+=2**(n-1-k)*comb(g,k)
    return h,q
def collisions():
    for n in range(2,10):
        seen={}
        for inn in product('AB',repeat=n-2):
            w='A'+''.join(inn)+'B'
            k=invariant(w)
            p=cycleprofile(w)
            if k in seen and p!=seen[k][1]:
                print(k,(w,p),seen[k])
                return w,seen[k][0]
            seen[k]=w,p
if __name__=='__main__':
    a,b=collisions()
    for s in range(9):
        u,v=(('A'*s)+w+('B'*s) for w in (a,b))
        print(s,invariant(u),cycleprofile(u),cycleprofile(v))
