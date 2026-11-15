from itertools import combinations,product

def profiles(n,code,R=3):
    V=1<<n; full=V-1
    O=[0]*n
    for j,(a,b) in enumerate(combinations(range(n),2)):
        if code>>j&1: O[a]|=1<<b
        else: O[b]|=1<<a
    h=[0]*V;h[0]=1;E=[[0]*n for i in range(V)]
    for m in range(1,V):
        s=m
        while s:
            t=s&-s;s^=t;v=t.bit_length()-1
            pred=(m^t)&~O[v]
            if m==t:E[m][v]=1
            while pred:
                i=pred&-pred;pred^=i;E[m][v]+=E[m^t][i.bit_length()-1]
        h[m]=sum(E[m])
    h[0]=0
    p=[h]
    for k in range(2,R+1):
        g=[0]*V
        for m in range(1,V):
            s=m
            while s:
                g[m]+=h[s]*p[-1][m^s]
                s=(s-1)&m
        p.append(g)
    return tuple(q[-1] for q in p)

def C3(blocks):
    lengths=[z[0] for z in blocks]; edge={}
    of=0
    for i,(n,c) in enumerate(blocks):
        for k,(a,b) in enumerate(combinations(range(n),2)):edge[a+of,b+of]=(c>>k)&1
        for a in range(of):
            for b in range(of,of+n):
                edge[a,b]=not (i==2 and a<lengths[0])
        of+=n
    return of,sum(edge[a,b]<<i for i,(a,b) in enumerate(combinations(range(of),2)))

def Q(p,k):
    from math import comb
    return sum(comb(k,i)*v for i,v in enumerate(p,1))

def eq(a,b):
    if a!=b: raise AssertionError((a,b))

if __name__=="__main__":
    for m,n in [(1,1),(2,2),(3,2),(3,3),(4,3)]:
        for a,c in [(0,0),(5,2),(29,10)]:
            p,q=profiles(m,a,R=4),profiles(n,c,R=4)
            u=profiles(*C3([(m,a),(n,c),(1,0)]),R=3)
            for k in (1,2,3):
                f=[Q(p,k+1)-Q(p,k),Q(p,k),Q(p,k)-Q(p,k-1)]
                g=[Q(q,k+1)-Q(q,k),Q(q,k),Q(q,k)-Q(q,k-1)]
                eq(k*sum(x*y for x,y in zip(f,g)),Q(u,k))
    print("45 mixed-slot checks passed")
    parts=[(5,76),(6,88),(5,12),(6,408)]
    a,b,c,d=[[Q(profiles(n,t),k) for k in (1,2,3)] for n,t in parts]
    eq((a,b,c,d),([15,120,525],[33,360,2025],[11,100,465],[45,432,2277]))
    def composite(A,B):
        return [a*b for a,b in zip(A,B)]
    def eq_future(T,h):
        a,b,c=T
        eq(h, (a+b)+2*(c+2*b-a))
    A,B=composite(a,b),composite(c,d)
    eq(A[:2],B[:2])
    eq((A[2],B[2]),(1063125,1058805))
    eq_future(A,2341755)
    eq_future(B,2333115)
    print("same-size two-slot collision and second-wrap separation passed")

    def sign(n,i,j):
        if n==1: return 1,0
        m=n//2
        if i<m and j<m: return sign(m,i,j)
        if i<m:
            s,a=sign(m,j-m,i); return s,a+m
        if j<m:
            s,a=sign(m,i-m,j); return s*(1 if j==0 else -1),a+m
        s,a=sign(m,j-m,i-m); return -s*(1 if j==m else -1),a
    eq([sign(16,i,j) for i in (1,10) for j in (4,15)],
       [(1,5),(-1,14),(-1,14),(1,5)])
    print("zero-divisor witness passed")

    def z(j):
        return (j**3+2*j)**2
    for k in range(1,81):
        L=9; M=1<<L
        if k%2:
            d=k if k%8 in (3,5) else k**3+2*k
            S={3*pow(k,2*L,M)*pow(d,2*t,M)%M for t in range(M//8)}
            target=set(range(3,M,8))
        elif k%4==2:
            a,b=k-1,k+1
            S={(pow(a,2*(L+t),M)*pow(z(a),u,M)
                 +pow(b,2*(L+t),M)*pow(z(b),u,M))%M
               for t in range(M//8) for u in range(M//8)}
            target=set(range(2,M,8))
        else:
            a,b=k-1,k+1
            S={(pow(a,2*L,M)*pow(z(a),t,M)
                 +pow(b,2*L,M)*pow(z(b),t,M))%M
               for t in range(M//16)}
            target=set(range(2,M,16))
        eq(S,target)
    print("all expected congruence sheets realized at 9 bits for k=1..80")
