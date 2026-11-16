"""Exact masks, independent literal head paths and small subset DP."""
from itertools import combinations

def K(d):
    return 2**(d-1)+(d==1)

def rev(x,n):
    return int(f'{x:0{n}b}'[::-1],2)

def J(s):
    return sum(K(j-i) for i,j in combinations(sorted(s),2))

def bits(x):
    return {i for i in range(x.bit_length()) if x>>i&1}

def pack(k,r,t,a,b):
    S={k-1-i for i in bits(b)}|{k+t+i for i in bits(a)}
    return k+r+t,S

def formula(k,r,t,a,b):
    P=(2**k-b)*(2**r+a)
    C=1+rev(b,k)-rev(a,r)-J(bits(a))-J(bits(b))
    return P,C,P*2**t+C

def head(N,marks):
    # Ascending tape including sentinel at N. Its inversions across
    # the apex are grouped by the old predecessor/successor of the head.
    Z=marks|{N}
    return 1+sum(K(i-j) for i in Z for j in range(i) if j not in Z)

def dp(N,marks):
    n=N+2;F=(1<<n)-1
    a=[((1<<(N+1))-1)^((1<<(i+1))-1) for i in range(N+1)]+[0]
    for i in marks|{N}:a[i]|=1<<(N+1)
    for i in range(N):
        if i not in marks:a[-1]|=1<<i
    T=[[0]*n for _ in range(F+1)]
    for i in range(n):T[1<<i][i]=1
    for s in range(1,F):
        for i,z in enumerate(T[s]):
            if not z:continue
            b=a[i] & ~s
            while b:
                j=b&-b;b^=j
                T[s|j][j.bit_length()-1]+=z
    return sum(T[F])

if __name__=='__main__':
    for r in range(5,11):
        D=1<<r
        for a in range(4,D//2-1,4):
            b=D-a;k=r+1
            C=formula(k,r,2,a,b)[1]
            P=(D+a)**2-1
            v=(a&-a).bit_length()-1
            if r>=v+3:
                assert C%4==(1-2*(a>>(v+2)&1))%4
            if a==3<<v:
                assert C==1
            for o,q,c in [(1,P,C),(-1,P,C+2),(0,P+1,C)]:
                h=formula(k,r,2,a+o,b+o)
                assert h==(q,c,4*q+c)
                N,V=pack(k,r,2,a+o,b+o)
                assert head(N,V)==h[-1]
                if r==5 and a==12:
                    assert dp(N,V)==h[-1]
                    print('radix4',o,h)
            N,V=pack(k,r,2,a-1,b-1)
            assert head(N,V|{k,k+1})==4*(P+1)+C
    for r in range(4,64):
        D=1<<r; k=r+1; t=r-1
        a=((D+1) if r%2 else (D-7))//3
        b=D-a
        P=(D+a)**2
        C=((6*r+7-2*D) if r%2 else (6*r+1-D))//9
        for o,q,c in [(-1,P-1,C),(1,P-1,C+D//2-2),(0,P,C)]:
            h=formula(k,r,t,a+o,b+o)
            assert h==(q,c,q*2**t+c)
            N,V=pack(k,r,t,a+o,b+o)
            assert head(N,V)==h[-1]
            if r==4:assert dp(N,V)==h[-1]
        if r%2 and r>=7:
            N,V=pack(k,r,t,a,b)
            for d in (1,3):
                j=(t-1-d)//2
                assert head(N,V|{k+j,k+j+d})==P*2**t+C-K(d)
    for r in range(4,10):
        D=1<<r; k=r+1; t=r-1
        a=((D+1) if r%2 else (D-7))//3
        b=D-a
        full=set(range(t))
        N,VH=pack(k,r,t,a,b)
        _,VL=pack(k,r,t,a-1,b-1)
        _,VT=pack(k,r,t,a+1,b+1)
        H=head(N,VH)
        palette=set()
        for x in range(1<<t):
            U=bits(x)
            W={t-1-j for j in full-U}
            lh=head(N,VL|{k+j for j in U})
            hh=head(N,VH|{k+j for j in W})
            assert lh==hh
            palette|={lh,hh,head(N,VT|{k+j for j in U})}
        if r>=7 and r%2:
            assert {H-4,H-2,H}<=palette and H-6 not in palette
    print('Whole spacer palette correspondence and immediate gap verified through r=9.')
    print('Radix-four coalescences and variable-height carries verified.')
