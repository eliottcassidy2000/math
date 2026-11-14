from collections import Counter
from itertools import permutations
def mul(p,q):
    r=Counter()
    for a,b in p.items():
        for c,d in q.items():r[a+c]+=b*d
    return {i:c for i,c in r.items() if c}
def det(M):
    n=len(M); acc=Counter()
    for q in permutations(range(n)):
        t={0:1-2*(sum(q[i]>q[j] for i in range(n) for j in range(i+1,n))%2)}
        for i in range(n):t=mul(t,M[i][q[i]])
        acc.update(t)
    return {i:c for i,c in acc.items() if c}
M=[[({0:1} if u==v else {})
    if (u<3)==(v<3) else ({u*(v-3):1} if u<3 else {-(u-3)*v:1})
    for v in range(6)] for u in range(6)]
assert det(M)=={-3:-4,3:-4}
norm=[]
for j in range(6):
    col=[]
    for i in (0,2,3,5):
        a=det([[M[u][v] for v in range(6) if v!=i] for u in range(6) if u!=j])
        col.append(sum(abs((x+2)%4-2) for x in a.values()))
    norm.append(sum(col))
assert norm==[20]*6
for r in (3,6,9,12,15,18,24,36):
    D={(0,j) for j in range(r) if j%3!=2}
    D|={(3+d,t) for d in range(3) for t in range(r) if t%3==2}
    def neighbors(v,k):
        if v<3:return { (3+d,(k-d*v)%r) for d in range(3)}
        return {(i,(k+(v-3)*i)%r) for i in range(3)}
    assert all((v,k) in D or neighbors(v,k)&D for v in range(6) for k in range(r))
    assert len(D)==5*r//3
    print(r, 'bounds', (8*r+4)//5, len(D))
print('Determinant and all six adjugate column norms verified')
