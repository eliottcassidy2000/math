from replay_intercyclic import rows
from itertools import permutations,product

def counts(w,z):
    R=rows(w)
    n=len(R);z+=2
    total=xy=xz=yz=0
    for p in permutations(range(n)):
        if not all(R[a]>>b&1 for a,b in zip(p,p[1:])):continue
        a,b,c=map(p.index,(0,1,z))
        total+=1;xy+=a<b;xz+=a<c;yz+=b<c
    return total,xy,xz,yz
for a,b in product((1,2,3),repeat=2):
    A,B=2**a,2**b;H=1+3*A*B
    assert counts('1'*a+'2'+'1'*b,a)==(H,H-A-B-1,2*A*B,A*B+1)
for a in (1,2,3):
    A=2**a;H=9*A-3
    assert counts('1'*a+'20',a)==(H,6*A-3,6*A-2,3*A)
print('Marked balance grids for both all-order families passed')
