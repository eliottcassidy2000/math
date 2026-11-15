from replay_recontact_energy import st, paths, ordinal, cap, reach
from math import isqrt

def p(n):
    v=0
    while n:
        q=n//2
        if n%2==0:v+=q+(q%2)
        n=q
    return v

def one(Y,B,h,e):
    out=[]
    if Y<0:return out
    q=2
    while (q*B)**2<h*Y:q*=2
    r=q.bit_length()-1
    for s in range(q):
        R=int(f'{s:0{r}b}'[::-1],2)
        A=q*B+h*(q-1-R)
        Z=Y-B*s-h*p(s)-e*(s%2)
        if Z<0:continue
        j=max(0,-((h-Z)//(A+h)))
        for u in range(j,Z//A+1):
            z=q*u+s
            if B*z+h*p(z)+e*(z%2)==Y:out.append(z)
    return out
def N(i,x,y):
    return 41+2*i+15*(p(x)+p(y))+63*x+42*y+120*x*y+9*(x%2)+6*(y%2)
def inverse(a):
    ret=set()
    for i in (0,1):
        for x in range(isqrt(a//120)+1):
            yl=a-41-2*i-15*p(x)-63*x-9*(x%2)
            for y in one(yl,42+120*x,15,6):ret.add((i,x,y))
        for y in range(isqrt(a//120)+1):
            xl=a-41-2*i-15*p(y)-42*y-6*(y%2)
            for x in one(xl,63+120*y,15,9):ret.add((i,x,y))
    return ret
def main():
    for r,n in [((16,1,3,7,14),25),((20,1,2,7,14),23)]:
        h,N0,s,d=st(r,22)
        vals=paths(r)[1];f=len(vals)-1
        Q=sum(x*vals[f^j] for j,x in enumerate(vals))
        if (h,Q,N0,s,d)!=(9,84,n,41,7):raise AssertionError((r,N0))
        T=ordinal(r,(0,))
        if st(T,22)!=(9,n+48,82,0) or not all(reach(cap(T,22),i) for i in range(len(T)+1)):raise AssertionError(T)
        print('one reversed edge profile and guard:',(h,Q,N0,s,d),n+48)
    L=2300
    expect={}
    for i in (0,1):
        for x in range(L//63+1):
            for y in range(L//42+1):
                a=N(i,x,y)
                if a<=L:expect.setdefault(a,set()).add((i,x,y))
    for a in range(1,L+1):
        q=inverse(a)
        if q!=expect.get(a,set()):raise AssertionError((a,q,expect.get(a)))
    for i,x,y in [(1,123,413),(0,1000,19),(0,0,3602),(1,3200,0),(0,78,79)]:
        if (i,x,y) not in inverse(N(i,x,y)):raise AssertionError((x,y))
    print('Exact inverse over both tails: all targets through',L,'and five larger controls')
if __name__=='__main__':main()
