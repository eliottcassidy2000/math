"""Busch's shared-vertex structural reduction, followed by literal path counts."""
from itertools import permutations, product

def H(n,E):
    out=[0]*n
    for a,b in E:
        if out[a]>>b&1 or out[b]>>a&1:
            raise ValueError((a,b))
        out[a]|=1<<b
    if sum(x.bit_count() for x in out)!=n*(n-1)//2:
        raise ValueError(out)
    return sum(all(out[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(n)))

def arcs(rows,labels):
    return [(a,b) for i,a in enumerate(labels) for j,b in enumerate(labels) if rows[i]>>j&1]

def seam(A,B,s):

    right=list(range(3,2+len(B)))
    labels=right.copy()
    labels.insert(s,0)
    e=arcs(A,(0,1,2))+arcs(B,labels)
    links=list(product((1,2),right))
    return {H(2+len(B),e+[(b,a) if sign else (a,b) for (a,b),sign in zip(links,ss)])
            for ss in product((0,1),repeat=len(links))}

C=(2,4,1); D=(6,8,10,1)
z=seam(C,C,0)|{9}
print('Order 5:',sorted(z))
for i in (0,1):
    print('Order 6 root orbit',i,':',sorted(seam(C,D,i)|{17}))
