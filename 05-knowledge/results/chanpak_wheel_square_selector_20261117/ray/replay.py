from math import isqrt

BANK={
  1:((56,41,43,16,6,24),12,23),
  7:((50,8,11,33,14,22),6,33),
  9:((42,48,3,6,45,12),40,37),
 13:((40,13,25,32,107,6,47),40,45),
 17:((114,64,59,51,66,18,44),66,45),
 19:((20,21,40,51,36,3),33,37),
 23:((42,48,3,6,45,12),34,37),
 25:((92,109,24,64,10,93,20),24,41),
 41:((118,60,16,53,96,4,46),64,41),
}
def check(value, expected):
    if value != expected: raise AssertionError((value,expected))

def check_graph(a):
    n=len(a);full=(1<<n)-1
    for i,x in enumerate(a):
        check(x & ~(full^(1<<i)),0)
        for j in range(i):
            check(((x>>j)&1)+((a[j]>>i)&1),1)

def strong(a):
    full=(1<<len(a))-1
    for i in range(len(a)):
        r=1<<i
        for _ in a:
            for j,x in enumerate(a):
                if r>>j&1:r|=x
        if r!=full:return False
    return True

def ends(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)];h=[1]+[0]*((1<<n)-1)
    for i in range(n):d[1<<i][i]=1
    for m in range(1,len(d)):
        h[m]=sum(d[m])
        for i,c in enumerate(d[m]):
            if not c:continue
            z=a[i]&~m
            while z:
                b=z&-z;z^=b
                d[m|b][b.bit_length()-1]+=c
    return h,d

def profile(a,q):
    h,d=ends(a);_,t=ends(tuple(sum(((a[j]>>i)&1)<<j for j in range(len(a))) for i in range(len(a))))
    e=[sum(c for j,c in enumerate(row) if q>>j&1) for row in d]
    s=[sum(c for j,c in enumerate(row) if q>>j&1) for row in t]
    f=len(h)-1
    return h[f],sum(h[f^m]*s[m] for m in range(f+1)),sum(e[f^m]*s[m] for m in range(f+1))

def physical(a,q,b):
    n=len(a);m=len(b);mask=((1<<m)-1)<<(n+1)
    return tuple(x|(int(not(q>>i&1))<<n)|mask for i,x in enumerate(a))+(q,)+tuple((x<<(n+1))|(1<<n) for x in b)

def prime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))

def main():
    for K,(a,q,J) in BANK.items():
        check_graph(a);check(strong(a),True)
        check(profile(a,q),(J,99,J-K))
        for b in [(0,), (2,0), (6,4,0), (2,4,1), (2,4,9,3)]:
            X=physical(a,q,b); check_graph(X);check(strong(X),True)
            vals=ends(b)[0];f=len(vals)-1
            check(ends(X)[0][-1],99*sum(vals[i]*vals[f^i] for i in range(f+1))+K*vals[-1])
        print(K,profile(a,q))
    for k in range(1,100):
        check([K for K in (13,19,17) if (99*(1<<k)+K)%7==0],[(13,19,17)[k%3]])
    check([k for k in range(1,17) if prime((99*(1<<k)+(13,19,17)[k%3])//7)],[1,2,4,6,8,9,13,14])
    check([k for k in range(1,17) if prime(99*(1<<k)+23) and prime(99*(1<<k)+25)],[2,4,6,12,14])
    first25=[p for p in range(2,100) if prime(p)]
    check(len(first25),25)
    n=99*(1<<16)+23
    check(n,1301*4987);check(any(n%p==0 for p in first25),False)
    for K,p in zip(BANK,(5,53,3,7,29,59,61,31,5)):
        check(p in first25,True);check(pow(2,158340,p),1);check((99+K)%p,0)
    print('first-25 composite',n)
    print('all nine rails composite when k is a positive multiple of 158340')
if __name__=='__main__':
    main()
