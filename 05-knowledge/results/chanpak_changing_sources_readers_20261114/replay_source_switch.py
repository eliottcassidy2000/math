from itertools import combinations, product
from functools import cache

def test(cond, desc):
    if not cond: raise AssertionError(desc)

def tables(A):
    n=len(A); mask=(1<<n)-1
    @cache
    def p(s,v,rev=False):
        if s==1<<v:return 1
        t=s^(1<<v)
        return sum(p(t,w,rev) for w in range(n) if t>>w&1
                   and (A[v]>>w&1 if rev else A[w]>>v&1))
    @cache
    def h(s):
        return sum(p(s,v) for v in range(n) if s>>v&1) if s else 1
    def pointed(r,starting):
        e=o=0
        for s in range(1,mask+1):
            if s>>r&1:
                z=p(s,r,starting)*h(mask^s)
                if s.bit_count()%2:o+=z
                else:e+=z
        return h(mask),e,o
    def cover(): return h(mask),sum(h(s)*h(mask^s) for s in range(mask+1))
    return mask,h,pointed,cover

def tt(n): return tuple((1<<n)-(1<<(j+1)) for j in range(n))
def K(A,B):
    n=len(A); m=len(B)
    return (((1<<(n+1))-2),)+tuple((a<<1)|(((1<<m)-1)<<(n+1)) for a in A)+tuple((b<<(n+1))|1 for b in B)
def seam(S,R,w):
    n=len(S); m=len(R)
    rows=[a|(((1<<m)-1)<<n) for a in S]+[b<<n for b in R]
    rows[0]^=1<<(n+w); rows[n+w]^=1
    return tuple(rows)
def small(n):
    ij=list(combinations(range(n),2))
    for bits in range(1<<len(ij)):
        A=[0]*n
        for k,(i,j) in enumerate(ij):
            if bits>>k&1:i,j=j,i
            A[i]|=1<<j
        yield tuple(A)
def decode(N,a,s,shifts):
    for z in shifts:
        k=1
        while ((1<<k)-1)*(2*s-a)+z<=N:
            u=(1<<k)-1
            if (N-z)%u==0:
                x=(N-z)//u+a
                if x%s==0:
                    y=x//s
                    if y>=2 and y&(y-1)==0:yield k,y.bit_length()-1,z
            k+=1

def main():
    pool=[A for n in range(1,4) for A in small(n)]
    profiles={A:tables(A)[3]() for A in pool}
    for A,B in product(pool,repeat=2):
        ha,qa=profiles[A]; hb,qb=profiles[B]; u=ha*hb; D=qa-ha
        val=tables(K(A,B))[2](0,True)
        test(val==(2*u+D*(qb-hb),D*qb//2,u+D*qb//2),('source',A,B,val))
        adj=tables(K(B,A))[2](0,True)
        det=(qa*hb-qb*ha)//2
        test(tuple(i-j for i,j in zip(val,adj))==(0,det,det),('swap',A,B))
    readers=[
        ((26,60,57,48,32,1),2,15,7,12),
        ((14,20,24,18,1),1,15,8,11),
        ((26,60,57,48,32,1),0,15,9,10),
        ((2,4,1),0,3,1,2),
        ((4,13,8,1),1,3,0,3)
    ]
    count=0
    for R,w,a,e,o in readers:
        test(tables(R)[2](w,False)==(a,e,o),('reader',R,w))
        s=a+e+o
        cases=[(tt(k),tt(m)) for k,m in product(range(1,4), repeat=2)]
        cases.extend([((2,4,1),(2,4,1)),((2,4,1),tt(2)),(tt(2),(2,4,1))])
        for A,B in cases:
            ha,qa=tables(A)[3](); hb,qb=tables(B)[3]()
            T=seam(K(A,B),R,w)
            mask,H,_,_=tables(T)
            want=(2*a+2*e)*ha*hb+(qa-ha)*(s*qb-a*hb)
            test(H(mask)==want,('seam',R,w,A,B,want))
            count+=1
    A=(16,17,3,23,4)
    B=(56,57,35,4,44,8)
    test(tables(A)[3]()==(5,64) and tables(B)[3]()==(19,256),'Bézout block profiles')
    for C,D,(R,w,*_),want in [(A,B,readers[3],88021),(B,A,readers[4],88023)]:
        T=seam(K(C,D),R,w)
        mask,H,_,_=tables(T)
        test(H(mask)==want,('Bezout literal glued count',want))
    triple={((1<<k)-1)*(34*(1<<m)-15)+z for k,m in product(range(1,12), repeat=2)
            for z in (44,46,48)}
    for N in range(1,4000):
        test(bool(list(decode(N,15,34,(44,46,48))))==(N in triple),('inverse',N))
    print(f'{len(pool)**2} literal sources and swaps, 5 literal readers, {count} glued Hamilton counts, literal Bézout pair 88021/88023, inverse through 3999')
if __name__=='__main__':main()
