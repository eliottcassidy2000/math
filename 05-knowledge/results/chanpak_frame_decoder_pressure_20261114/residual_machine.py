"""Self-contained exact Hankel certificate for ordinary two-hub appends."""
from fractions import Fraction as Q
from collections import deque
from counter36 import AUTO, H


def scalar(a,b):
    return sum(x*y for x,y in zip(a,b))


class Basis:
    def __init__(self):
        self.pivots={}
    def add(self,v):
        u=list(map(Q,v))
        for i,p in self.pivots.items():
            if u[i]:
                a=u[i]
                for k in range(i,len(u)):
                    u[k]-=a*p[k]
        for i,a in enumerate(u):
            if a:
                self.pivots[i]=[x/a for x in u]
                return True
        return False


def rank(matrix):
    basis=Basis()
    return sum(basis.add(r) for r in matrix)


def closure(v,step,append):
    basis=Basis()
    out=[];todo=deque([(v,'')])
    while todo:
        u,s=todo.popleft()
        if basis.add(u):
            out.append((u,s))
            for d in range(4):
                todo.append((step(u,d),append(s,str(d))))
    return out


def solve(A,B):
    n=len(A); rows=[list(map(Q,a+b)) for a,b in zip(A,B)]
    for i in range(n):
        j=next(j for j in range(i,n) if rows[j][i])
        rows[i],rows[j]=rows[j],rows[i]
        a=rows[i][i]
        rows[i]=[x/a for x in rows[i]]
        for j in range(n):
            if j!=i and rows[j][i]:
                a=rows[j][i]
                rows[j]=[x-a*y for x,y in zip(rows[j],rows[i])]
    return [row[n:] for row in rows]


A,B=AUTO
n=len(A[0][0])
N=n+len(B[0][0])
T=[list(A[0][d])+[tuple(j+n for j in v) for v in B[0][d]] for d in range(4)]
f=tuple(int(i in A[1] or i-n in B[1]) for i in range(N))
e=tuple(int(i in (0,n)) for i in range(N))

def onward(v,d):
    u=[0]*N
    for i,a in enumerate(v):
        for j in T[d][i]:
            u[j]+=a
    return tuple(u)

def backward(v,d):
    return tuple(sum(v[j] for j in choices) for choices in T[d])

PROBES=('', '0','1','2','00','10','20','30','01','11','21','31',
        '000','100','200','010','210')
PREFIXES=('', '1','2','3','10','11','12','13','20','21','22','23',
          '102','103','120','122','123')

def profile(w):
    return tuple(H(w+z) for z in PROBES)

def advance(v,matrix):
    vals=tuple(sum(x*y for x,y in zip(v,col)) for col in zip(*matrix))
    assert all(x.denominator==1 for x in vals)
    return tuple(x.numerator for x in vals)

def certify():
    reach=closure(e,onward,lambda a,b:a+b)
    obs=closure(f,backward,lambda a,b:b+a)
    assert [N,len(reach),len(obs)]==[36,28,20]
    assert rank([[scalar(a,b) for b,_ in obs] for a,_ in reach])==17
    M=[list(profile(w)) for w in PREFIXES]
    assert rank(M)==17
    D=[]
    for d in '0123':
        mat=solve(M,[list(profile(w+d)) for w in PREFIXES])
        assert [17-rank([[a-(i==j)*ev for j,a in enumerate(row)]
                         for i,row in enumerate(mat)]) for ev in range(4)]==[4,8,4,1]
        for a,w in reach:
            assert advance(profile(w),mat)==profile(w+d)
        D.append(mat)
    return D


if __name__=='__main__':
    certify()
    assert H('10011')==H('1121')==25 and H('100110')==63 and H('11210')==79
    print('36 states; reachable 28; observable 20; Hankel rank 17; four letter maps all rank 13')
    print('same 25, distinct continuations:', H('100110'), H('11210'))
