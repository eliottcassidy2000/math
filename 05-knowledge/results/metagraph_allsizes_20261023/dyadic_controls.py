"""Independent bounded lift enumeration for higher dyadic energy forms."""
from itertools import product
import json

CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise ValueError(m)

def product_matrix(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def transpose(A):return list(map(list,zip(*A)))
def conjugate(A,U):return product_matrix(transpose(U),product_matrix(A,U))
def image(A,x):return [sum(a*b for a,b in zip(row,x)) for row in A]
def rank(rows):
    piv={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p in piv:x^=piv[p]
            else:piv[p]=x;break
    return len(piv)

def level(A,j):
    n=len(A);mod=1<<j;bank={}
    for x in product(range(mod),repeat=n):
        if all(v%mod==0 for v in image(A,x)):
            key=sum((v%2)<<i for i,v in enumerate(x));bank.setdefault(key,[]).append(x)
    basis=[]
    for x in sorted(bank):
        if rank(basis+[x])>len(basis):basis.append(x)
    need(len(bank)==1<<len(basis),'full binary subgroup')
    rows=[]
    for u in basis:
        row=0;x=bank[u][0]
        for t,v in enumerate(basis):
            y=bank[v][0];Ay=image(A,y);value=sum(a*b for a,b in zip(x,Ay))
            need(value%mod==0,'energy integral')
            bit=(value//mod)%2;row|=bit<<t
            for other in bank[u]:
                need((sum(a*b for a,b in zip(other,Ay))//mod)%2==bit,'all bounded lift choices agree')
        rows.append(row)
    radical=set()
    for coefficients in range(1<<len(basis)):
        if all((row&coefficients).bit_count()%2==0 for row in rows):
            u=0
            for i,b in enumerate(basis):
                if coefficients>>i&1:u^=b
            radical.add(u)
    return set(bank),radical,rank(rows)

def run():
    U=[[1,1,0],[0,1,1],[1,1,1]] # determinant exactly 1
    cases=[]
    for diag in [(2,4,8),(1,2,16)]:
        D=[[diag[i] if i==j else 0 for j in range(3)] for i in range(3)]
        cases += [(D,diag),(conjugate(D,U),diag)]
    for e in (1,2,3):
        cases.append(([[0,1<<e],[1<<e,0]],(1<<e,1<<e)))
    answers=[]
    for A,smith in cases:
        previous=None;rows=[]
        for j in range(1,5):
            V,R,r=level(A,j)
            need(len(V)==1<<sum(s%(1<<j)==0 for s in smith),'Smith divisibility dimension')
            need(len(R)==1<<sum(s%(1<<(j+1))==0 for s in smith),'radical next divisibility')
            if previous is not None:need(previous==V,'radical equals next liftable space')
            previous=R;rows.append([j,len(V),len(R),r])
        answers.append({'A':A,'known_smith':smith,'levels':rows})
    print(json.dumps({'checks':CHECKS,'cases':answers},separators=(',',':')))

if __name__=='__main__':run()
