"""Construct arbitrary rectangular certificates, and independently check the small path counts."""
from itertools import permutations
from math import isqrt

def need(c, *w):
    if not c: raise AssertionError(w)

def drive(a,b):
    def g(A,B):
        return A*B-min(A,B) if min(A,B)>=2 else A*B
    k=next(k for k in range(min(a,b)+1) if g(max(a-k,0),max(b-k,0))<=k)
    X=list(range(1,a+1));Y=list(range(a+1,a+b+1));v=0
    def arc(i,j):
        if i==j:return False
        if not i:return j in X
        if not j:return i in Y
        return i<j
    if k==min(a,b):
        if b<=a:
            order=[v]+X+Y
            paths=[Y[:i+1]+[v]+X+Y[i+1:] for i in range(b)]
        else:
            order=X+Y+[v]
            paths=[X[:i]+Y+[v]+X[i:] for i in range(a)]
    else:
        A=X[k:];D=Y[:-k]
        order=X[:k]+D+[v]+A+Y[-k:]
        cells={(x,y) for x in A for y in D}
        groups=[]
        if min(len(A),len(D))>=2:
            groups=[[(A[-1],y),(x,D[0])] for x,y in zip(A[:-1],D[1:])]
            groups.append([(A[-1],D[0]),(A[0],D[1])])
        for group in groups:
            for e in group:cells.remove(e)
        groups.extend([[e] for e in sorted(cells)])
        groups.extend([[] for i in range(k-len(groups))])
        paths=[]
        for es,u,z in zip(groups,X[:k],Y[-k:]):
            if not es:
                preX=[x for x in X if x<u]; preY=[y for y in Y if y<=z]
            elif len(es)==1:
                x,y=es[0]
                preX=[i for i in X if i<u]+[x];preY=[y,z]
            else:
                (x,y),(xx,yy)=es
                preX=[i for i in X if i not in [u,xx] and (i==x or i<u or i>xx)]
                preY=[j for j in Y if j not in [yy] and (j in [y,z] or j<yy)]
            p=preX+preY+[v]+[i for i in X if i not in preX]+[j for j in Y if j not in preY]
            paths.append(p)
    need(len(paths)==k,a,b,k,len(paths))
    seen=set()
    for p in paths:
        need(sorted(p)==list(range(a+b+1)),a,b,p)
        for i,j in zip(p,p[1:]):
            need(arc(i,j),a,b,p,i,j)
            seen.add((i,j))
    need(all((j,i) in seen for t,i in enumerate(order) for j in order[t+1:] if arc(j,i)),a,b,k)
    return arc,k

def main():
    for a in range(1,43):
        for b in range(1,43):
            arc,k=drive(a,b)
            if a==b:need(k==max(1,a-isqrt(a)),a,k)
            if a+b<=8:
                H=sum(all(arc(i,j) for i,j in zip(p,p[1:])) for p in permutations(range(a+b+1)))
                need(H==(2**a-1)*(2**b-1)+2,a,b,H)
    print('1764 literal covers; small independent permutation counts: PASS')

if __name__=='__main__':main()
