"""Emit rooted chart in the wedge screen input format."""
import sys
n=int(sys.argv[1]); assert n in [5,10,50]
A=[set() for _ in range(n)]
def add(u,v):A[u].add(v);A[v].add(u)
if n==5:
    for x in range(5):add(x,(x+1)%5)
else:
    q=n//10
    for i in range(q):
        for x in range(5):
            add(5*i+x,5*i+(x+1)%5)
            add(5*q+5*i+x,5*q+5*i+(x+2)%5)
            for j in range(q):add(5*i+x,5*q+5*j+(x+i*j)%5)
order=[0]+sorted(A[0])
for h in sorted(A[0]):order+=sorted(A[h]-{0})
assert len(set(order))==len(order)==n
I={a:i for i,a in enumerate(order)}
print(n)
for v in order:print(len(A[v]),*(sorted(I[w] for w in A[v])))
