from fractions import Fraction
from math import comb

def edge_orbits(n):
    edges=[(i,j) for i in range(n) for j in range(i+2,n)]
    left=set(edges); orbits=[]
    while left:
        e=min(left); r=(n-1-e[1],n-1-e[0]); o=sorted({e,r}); orbits.append(o);left.difference_update(o)
    return edges,orbits

def get_lap(n,orbits,mask):
    L=[[0]*n for _ in range(n)]
    for i in range(n-1): L[i][i]+=1;L[i][i+1]-=1
    for k,o in enumerate(orbits):
        for i,j in o:
            if mask>>k&1: a,b=i,j
            else:a,b=j,i
            L[a][a]+=1;L[a][b]-=1
    return L

def forest_coeff(L):
    n=len(L); P=[row[:] for row in L]; traces=[0]
    for k in range(1,n):
        traces.append(sum(P[i][i] for i in range(n)))
        if k+1<n:
            P=[[sum(P[i][k]*L[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    c=[1]
    for j in range(1,n):
        x=sum((-1)**(k-1)*c[j-k]*traces[k] for k in range(1,j+1));assert x%j==0;c.append(x//j)
    return c

if __name__=='__main__':
    for n in range(2,9):
        edges,orbs=edge_orbits(n);B=2**len(orbs);s=[0]*n
        for mask in range(B):
            f=forest_coeff(get_lap(n,orbs,mask));s=[x+y for x,y in zip(s,f)]
        all_expected=[Fraction(0)]*n
        for l in range((n-1)//2+1):
            for j in range(n-2*l):
                all_expected[2*l+j]+=Fraction(comb(n-1-l,l)*comb(n-1-2*l,j)*n**j,4**l*2**j)
        b=[Fraction(x,B) for x in s]
        print(n,'B',B,'blue',list(map(str,b)),'difference',list(map(str,[x-y for x,y in zip(b,all_expected)])),flush=True)
