from itertools import product

def counts(rows):
    n=len(rows); f=(1<<n)-1
    dp=[[0]*n for _ in range(f+1)]
    h=[1]+[0]*f
    for m in range(1,f+1):
        for i in range(n):
            if not(m>>i&1):continue
            s=m^(1<<i)
            dp[m][i]=1 if not s else sum(dp[s][j] for j in range(n) if s&rows[i]&(1<<j))
        h[m]=sum(dp[m])
    return h

def covers(h,k):
    f=len(h)-1;x=[1]+[0]*f
    for _ in range(k):
        y=[0]*(f+1)
        for m in range(f+1):
            s=m
            while True:
                y[m]+=x[m^s]*h[s]
                if not s:break
                s=(s-1)&m
        x=y
    return x[f]

def wrap(r):
    n=len(r)
    return tuple(x|(3<<n) for x in r)+(1<<(n+1)|3<<(n+2),3<<(n+2),((1<<n)-1)|1<<(n+3),(1<<n)-1)

def feasible(A,S,D):
    if A<0 or S<0 or A>S or (A-S)%3:return False
    for _ in range(D-1):
        S,d=divmod(S,4)
        A-=d
    return A>=S

def strong(rows):
    for i in range(len(rows)):
        reach=1<<i
        for _ in rows:
            for j,r in enumerate(rows):
                if reach>>j&1:reach|=r
        if reach!=2**len(rows)-1:return False
    return True

def main():
    for r,q3 in [((92,141,32,228,206,19,38,101),45225),((52,229,208,231,10,148,49,81),44505)]:
        h=counts(r); got=h[-1],covers(h,2),covers(h,3)
        w=counts(wrap(r))[-1]
        if not strong(r) or got!=(345,5640,q3) or w!=4*q3+5*5640-11*345:raise AssertionError((r,got,w))
        print('profile and literal cyclic wrap',got,w)
    for D in range(1,6):
        Smax=100
        seen={}
        def visit(d,s,a):
            if d==D:
                seen.setdefault(s,set()).add(a);return
            for i in range((Smax-s)//4**d+1):visit(d+1,s+i*4**d,a+i)
        visit(0,0,0)
        for s,a in product(range(Smax+1),repeat=2):
            if feasible(a,s,D)!=(a in seen.get(s,())):raise AssertionError((D,s,a))
    print('bounded-digit cones checked through five coordinates and sum 100')
    for n in range(5,16):
        D=(n-1)//2;r=n-2*D
        base=[n-2*d for d in range(1,D+1)]
        t=sum(b*4**(D-j-1) for j,b in enumerate(base))
        for h in range(2*sum(base)+3,500,2):
            A=(h-1)//2-sum(base)
            for q in range(max(2**n,(n+1)*h),2**n*(h+3)//4+1):
                k=q-2**n
                if k%(1<<(r+1)) or (q-(-1)**n*h)%3:continue
                S=k//(1<<(r+1))-t
                if A>=0 and S>=A and (S-A)%3==0 and not feasible(A,S,D):
                    print('strong profile passing elementary floor parity and congruence, rejected by digit mass',n,h,q,A,S)
                    return
    raise AssertionError('no hostile')
if __name__=='__main__':main()
