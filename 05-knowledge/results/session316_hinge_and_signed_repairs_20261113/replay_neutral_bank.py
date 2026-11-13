from functools import cache
from itertools import permutations, product

R=(24,17,3,6,12)
def check(cond):
    if not cond: raise AssertionError
def rev(A):
    return tuple(sum(1<<j for j,a in enumerate(A) if a>>i&1) for i in range(len(A)))
def compiler(A,U):
    n=len(A); F=(1<<n)-1
    @cache
    def paths(S,i,j):
        if not S>>i&1 or not S>>j&1:return 0
        if S==1<<i:return int(i==j)
        if i==j:return 0
        return sum(paths(S^(1<<j),i,k) for k in range(n) if A[k]>>j&1)
    @cache
    def h(S):
        return sum(paths(S,i,j) for i in range(n) for j in range(n)) if S else 1
    @cache
    def q(S,t):
        if not t:return int(not S)
        return sum(h(B)*q(S^B,t-1) for B in range(F+1) if B&S==B)
    @cache
    def s(S):return sum(paths(S,i,j) for i in range(n) if U>>i&1 for j in range(n))
    @cache
    def e(S):return sum(paths(S,i,j) for j in range(n) if not U>>j&1 for i in range(n))
    def g(t):
        return (sum(s(B)*q(F^B,t) for B in range(F+1)),
                sum(e(B)*q(F^B,t) for B in range(F+1)),
                sum(s(B)*e(C)*q(F^B^C,t) for B in range(1,F+1) for C in range(1,F+1) if not B&C))
    return q,g
def extend(A,y):
    return tuple(a|(1<<len(A) if i==y else 0) for i,a in enumerate(A))+(((1<<len(A))-1)^(1<<y),)
def profile(A,U):
    q,g=compiler(A,U);F=(1<<len(A))-1
    s,e,b=g(0);S,E,_=g(1)
    return q(F,1),q(F,2),s,S,e,E,b
def fuse(X):
    from math import prod
    n=len(X)
    h,q,s,S,e,E,b=map(list,zip(*X))
    def outside(J):return prod(x for i,x in enumerate(h) if i not in J)
    return (s[0]*outside({0})+e[-1]*outside({n-1})+
        sum(b[i]*outside({i}) for i in range(n))+
        sum(e[i]*s[i+1]*outside({i,i+1}) for i in range(n-1))+
        sum(S[j]*E[i]*prod(q[j+1:i])*outside(set(range(j,i+1))) for j in range(n) for i in range(j+1,n)))
def forms(word,k=None,v=0):
    blocks=[(R,u) for u in word]+[(rev(R),31^u) for u in word[::-1]]
    if k is not None:blocks[k]=(extend(R,v),blocks[k][1])
    return [profile(A,U) for A,U in blocks]
for U in (3,5):
    q,g=compiler(R,U)
    for t in range(6):
        P=t**4+2*t**3+8*t*t+7*t+3
        check(g(t)==(2*P,3*P,2*(t+1)*(3*t*t+6*t+13+(U==5))))
        check(q(31,t)==t*(t**4+10*t*t+4))
    for v in range(5):
        reach=U&~(1<<v)
        for _ in range(5):
            for j in range(5):
                if reach>>j&1:reach|=R[j]&~(1<<v)
        check(reach==31^(1<<v))
for v,delta in ((0,(1,1,-1,-5,0)),(1,(-2,-13,-3,-22,-15))):
    a=profile(extend(R,v),3);b=profile(extend(R,v),5)
    check(a[:2]==b[:2]==(33,360))
    check(tuple(y-x for x,y in zip(a[2:],b[2:]))==delta)

for d in range(2,5):
    n=2*d
    for k in range(d):
        w=[3]*d;j=next(i for i in range(d) if i!=k);w[j]=5
        z=w.copy();z[k]=5;z[j]=3
        check(fuse(forms(w))==fuse(forms(z)))
        actual=fuse(forms(z,k))-fuse(forms(w,k))
        c=15 if not k else 9
        c-=(6 if k+1<d else 9)+66
        S=sum((63 if i<d else 42)*8**(i-k-1) for i in range(k+1,n))
        F=210*sum(8**r for r in range(k))
        check(actual==(c+S-F)*15**(n-2)>0)
        check(fuse(forms(z,k,1))<fuse(forms(w,k,1)))
check(fuse(forms([3,5]))==33787125)
check(fuse(forms([3,5],0))==74776365)
check(fuse(forms([5,3],0))==75458115)
print('slot polynomial, strong deletion, and signed probe replays passed')
