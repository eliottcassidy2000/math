"""Finite end of an all-order proof; see RECTANGLE.md."""
from itertools import permutations, product

def h(x):
    ans=1
    while x:
        x,t=divmod(x,2)
        if t==0:
            ans+=2*((x+1)//2)
    return ans

def J(a,b):return h(a)+h(b)-1+a*b+(a&1)*(b&1)

checks=0
for a,b,A,B in product(range(18),range(18),range(1,10),range(1,10)):
    if A*B>9:continue
    r={J(a,b),J(a+A,b),J(a,b+B),J(a+A,b+B)}
    if len(r)==4 and max(r)-min(r)<=6:
        raise AssertionError((a,b,A,B,r))
    checks+=1

def brute(w):
    n=len(w)
    def arc(i,j):
        if i==n:return w[j]=='1'
        if j==n:return w[i]=='0'
        return i<j
    return sum(all(arc(i,j) for i,j in zip(p,p[1:]))
               for p in permutations(range(n+1)))

left=(0,1)
right=(5,6)
for a,b in product(left,right):
    V=''.join(str(1-int(x)) for x in f'{b:03b}'[::-1])
    if brute(str(a)+V)!=J(a,b):
        raise AssertionError('literal rectangle mismatch')
witness=tuple(J(a,b) for b in right for a in left)
if witness!=(3,9,5,11):
    raise AssertionError(witness)
print(f'{checks} exact residual rectangles checked; sharp literal span: {witness}')
