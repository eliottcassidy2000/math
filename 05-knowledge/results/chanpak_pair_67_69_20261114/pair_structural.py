#!/usr/bin/env python3
from itertools import combinations, product, permutations
A=(2,12,17,37,75,23,47)
B=(2,12,17,37,43,71,31)
C3=((4,1,2),1,2,(1,2,0))
C5=((6,12,24,17,3),1,16,(4,0,1,2,3))
C7=((14,28,104,112,37,67,19),9,68,(2,3,5,6,4,0,1))

def tt(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def ends(r):
    n=len(r);F=(1<<n)-1
    e=[[0]*n for _ in range(F+1)]
    for i in range(n):e[1<<i][i]=1
    for m,z in enumerate(e):
        for i,t in enumerate(z):
            if t:
                b=r[i]&~m
                while b:
                    v=b&-b;b^=v;e[m|v][v.bit_length()-1]+=t
    return e

def H(r):return sum(ends(r)[-1]) if r else 1

def form(C,U,V,w,g,hi):
    """C exterior, W source tape, G sink tape; w,g are pairs (rows, contact mask)."""
    w,b=w;g,c=g
    N=len(w)+len(C)+len(g)+2
    z=[0]*N;s=N-2;t=N-1;K=len(w);P=K+len(C)
    for off,r in [(0,w),(K,C),(P,g)]:
        for i,a in enumerate(r):
            z[i+off]|=a<<off
            z[i+off]|=((1<<s)-1)^((1<<(off+len(r)))-1)
    for i in range(s):
        if i<K or i>=P:
            x=(b>>i&1) if i<K else c>>(i-P)&1
            if x:z[i]|=(1<<s)|(1<<t)
            else:z[s]|=1<<i;z[t]|=1<<i
        else:
            if U>>(i-K)&1:z[s]|=1<<i
            else:z[i]|=1<<s
            if V>>(i-K)&1:z[t]|=1<<i
            else:z[i]|=1<<t
    if hi:z[s]|=1<<t
    else:z[t]|=1<<s
    return tuple(z)

def tape(word):return tt(len(word)),sum(int(x)<<i for i,x in enumerate(word))
def form3(w,g,hi):
    C,U,V,_=C3
    return form(C,U,V,tape(w),tape(g),hi)

def verify(C,U,V,p):
    n=len(C);F=(1<<n)-1
    assert all((C[i]>>j&1)==(C[p[i]]>>p[j]&1) for i in range(n) for j in range(n))
    assert sum(1<<p[i] for i in range(n) if U>>i&1)==V
    e=ends(C);s=ends(tuple(F^r^(1<<i) for i,r in enumerate(C)))
    def D(U,V):
        return sum(sum(e[m][i] for i in range(n) if U>>i&1)*
                    sum(s[F^m][i] for i in range(n) if V>>i&1) for m in range(1,F))
    assert D(V,U)-D(U,V)==2
    for w,g in [('', ''),('', '01'),('01','101'),('0011','0101')]:
        a,b=H(form(C,U,V,tape(w),tape(g),False)),H(form(C,U,V,tape(w),tape(g),True))
        assert b-a==2
        print('exterior',n,'W',w or '-','G',g or '-','pair',a,b)
def main():
    for args in (C3,C5,C7):verify(*args)
    assert H(form3('01','101',False))==1489
    assert H(form3('01','101',True))==1491
    # Non-transitive replacements reveal the product factor in the same proof.
    W=((2,4,1),5);G=((2,4,1),2)
    for C,U,V,_ in (C3,C5,C7):
        assert H(form(C,U,V,W,G,True))-H(form(C,U,V,W,G,False))==18
    # The reflected ordinary pair and both prime modular quotients.
    assert (H(A),H(B))==(67,69)
    def module(T,S):
        return all(T[i]&S in (0,S) for i in range(len(T)) if not S>>i&1)
    assert all(not module(T,S) for T in (A,B) for S in range(1,127) if S.bit_count()>1)
    def reversors(T):
      return [p for p in permutations(range(7)) if
              all((T[i]>>j&1)==(T[p[j]]>>p[i]&1) for i in range(7) for j in range(i))]
    assert reversors(A)==[] and reversors(B)==[tuple(reversed(range(7)))]
    def blow(r,i,n):
        out=[j for j in range(len(r)) for _ in range(n if j==i else 1)]
        return tuple(sum(1<<b for b,j in enumerate(out)
            if (u!=j and r[u]>>j&1) or (u==j and a<b)) for a,u in enumerate(out))
    for r in range(1,6):
        assert H(blow(A,6,r))==34*(1<<r)-1 and H(blow(B,6,r))==32*(1<<r)+5
    print('PASS: exterior skew, two painted banks, five clocks, nontransitive controls and converse.')
if __name__=='__main__':main()
