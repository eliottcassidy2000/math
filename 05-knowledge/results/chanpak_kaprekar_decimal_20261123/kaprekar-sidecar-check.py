"""Exact standard-library replays for Kaprekar borrows, orbit and physical seeds."""
from collections import Counter
from itertools import combinations_with_replacement, product
from math import comb, gcd

def digits(x,B,n):
    out=[]
    for i in range(n):x,a=divmod(x,B);out.append(a)
    assert x==0
    return tuple(reversed(out))

def value(d,B):
    v=0
    for a in d:v=B*v+a
    return v

def gaps(d):
    z=sorted(d,reverse=True)
    return tuple(z[i]-z[-1-i] for i in range(len(z)//2))

def borrow(c,B,n):
    k=sum(x>0 for x in c)
    if k==0:return (0,)*n
    return c[:k-1]+(c[k-1]-1,)+(B-1,)*(n-2*k)+tuple(B-1-c[i] for i in range(k-1,0,-1))+(B-c[0],)

def fibre(c,B,n):
    if not c:return B
    z=(B-c[0])*(c[-1]+1)**(n%2)
    for a,b in zip(c,c[1:]):z*=a-b+1
    return z

def inverse_binary(y):
    if y<=0 or y%2==0:return None
    n=y.bit_length()
    if y+1==1<<n:return (1,n)
    s=(1<<n)+1-y
    bits=[i for i in range(s.bit_length()) if s>>i&1]
    if len(bits)==1:a=b=bits[0]-1
    elif len(bits)==2:a,b=bits
    else:return None
    if a>=2 and a+b==n:return a,b

def endpoints(rows):
    n=len(rows)
    assert all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i in range(n) for j in range(i))
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for s,row in enumerate(d):
        for i,v in enumerate(row):
            if not v:continue
            z=rows[i]&~s
            while z:
                b=z&-z;z-=b;d[s|b][b.bit_length()-1]+=v
    return d

def signed_shell(s):
    A=(8,1,3,166,15,23,63,119)
    n=len(A); F=(1<<n)-1; I=80
    D=endpoints(A); R=endpoints([r^F^(1<<i) for i,r in enumerate(A)])
    G=L=0
    for a in range(1,F+1):
        b=F^a
        z=sum(x for i,x in enumerate(D[a]) if I>>i&1)
        G+=z*(sum(R[b]) if b else 1)
        L+=z*sum(x for i,x in enumerate(R[b]) if I>>i&1)
    assert (G,sum(D[-1])-L)==(224,-7)
    f=((1<<s)-1)<<(n+1)
    return [r|((I>>i&1)<<n) for i,r in enumerate(A)]+[(F^I)|f]+[
        F|(f&~((1<<(i+1))-1)) for i in range(n+1,n+1+s)]

if __name__=="__main__":
    for B,n in product(range(2,9),range(1,8)):
        tally=Counter()
        for z in combinations_with_replacement(range(B),n):
            c=gaps(z)
            x=value(z[::-1],B)-value(z,B)
            assert digits(x,B,n)==borrow(c,B,n)
            assert x%(B-1)==0
            assert x%2==(c[0]%2 if B%2==0 and c else 0)
            tally[x]+=1
        assert len(tally)==comb(n//2+B-1,B-1)
        for z in combinations_with_replacement(range(B),n):
            c=gaps(z); x=value(borrow(c,B,n),B)
            assert tally[x]==fibre(c,B,n)
    seen={}
    for n in range(2,129):
        for a in range(1,n//2+1):
            b=n-a; h=((1<<a)-1)*((1<<b)-1)
            assert h not in seen and inverse_binary(h)==(a,b)
            assert digits(h,2,n)==(1,)*(a-1)+(0,)+(1,)*(b-a)+(0,)*(a-1)+(1,)
            assert sum(digits(h,2,n))==b
            seen[h]=n
    for L in range(1,65):
        assert sum(x<1<<L for x in seen)==L*L//4+1
    for B in range(2,102,2):
        for a in range(1,B):
            z=abs(2*a-(B+1)); v=z; period=0
            while True:
                v=abs(2*v-(B+1));period+=1
                if v==z:break
            q=(B+1)//gcd(a,B+1)
            assert pow(2,period,q) in (1,q-1)
            assert all(pow(2,j,q) not in (1,q-1) for j in range(1,period))
    for H,A in {35:(24,117,121,114,96,65,1),
                49:(106,52,121,114,97,64,2),
                63:(62,28,24,48,32,70,159,63),
                217:(70,60,24,17,33,77,158,63)}.items():
        assert sum(endpoints(A)[-1])==H
    for s in range(1,6):assert sum(endpoints(signed_shell(s))[-1])==7*((1<<(s+5))-1)
    print("Exact borrow, fibre, orbit, binary inverse and physical-seed controls passed.")
