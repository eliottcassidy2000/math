from itertools import product
from functools import cache

def rev(x,k):
    return int(format(x,f'0{k}b')[::-1],2)
def rc(x,k):
    return (1<<k)-1-rev(x,k)
def common(x,k):
    return tuple(3*int(z) for z in format(x,f'0{k}b'))
def values(w):
    a,b=zip(*(divmod(z,2) for z in w));n=len(w)
    la=sum(x<<i for i,x in enumerate(a[::-1]))
    lb=sum(x<<i for i,x in enumerate(b[::-1]))
    ra=sum((1-x)<<i for i,x in enumerate(a))
    rb=sum((1-x)<<i for i,x in enumerate(b))
    def F(j,i):
        return (1-a[i])*b[j]-(1-b[i])*a[j]
    d=b[0]-a[0]+b[-1]-a[-1]+sum((1<<(i-j-1))*F(j,i) for i in range(n) for j in range(i))+sum(F(j+1,j) for j in range(n-1))
    return d,la*rb-lb*ra

def literal(w,ab):
    n=len(w)
    rows=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]+[0,0]
    for i,z in enumerate(w):
        for a,v in zip(divmod(z,2),[n,n+1]):
            if a:rows[v]|=1<<i
            else:rows[i]|=1<<v
    rows[n+int(not ab)]|=1<<(n+int(ab))
    @cache
    def f(M,i):
        if M==1<<i:return 1
        S=M^(1<<i)
        return sum(f(S,j) for j in range(n+2) if (S&rows[i])&(1<<j))
    return sum(f((1<<(n+2))-1,i) for i in range(n+2))

for j in range(1,9):
    k=2*j+1
    g=((1<<(2*j))-1)//3
    seen=set()
    for e in product((0,1),repeat=j):
        X=(1<<(k-1))+sum(ei<<(2*i) for i,ei in enumerate(e))
        Y=X+g
        E=sum(ei<<(2*(j-1-i)) for i,ei in enumerate(e))
        assert Y<1<<k and rc(Y,k)-rc(X,k)==-4*g+6*E
        w=common(X,k)+(0,1,2,3)+common(rc(Y,k),k)
        delta,D=values(w)
        assert (delta,D)==(2*g,(1<<(k+1))*(-16*(1<<k)*g-4*g+6*E))
        assert D<0
        if j==1:
            assert literal(w,True)-literal(w,False)==delta
        seen.add(D)
    assert len(seen)==1<<j
for k,m in product(range(2,10),range(1,7)):
    s=common((1<<(k-1))-1,k)
    p=s+(0,1,2,3)+s[::-1]
    u=(3,)+(0,)*(m-1)
    w=u+p+tuple(3-z for z in u[::-1])
    assert values(w)==(2,-(1<<(k+1+2*m))*(15*(1<<k)+3))
print('510 distinct carry menus checked, plus 48 fixed-gap wrappers and literal small controls.')

def fib(n):
    a,b=0,1
    for _ in range(n):a,b=b,a+b
    return a

def menu_size(k,g):
    a,b=1,0
    for i in range(k-1):
        if g>>i&1:b+=a
        else:a+=b
    return a

for k in range(2,10):
    best=0
    for g in range(1<<(k-1)):
        found=set()
        for x in range(1<<(k-1),(1<<k)-g):
            y=x+g
            c=sum(int((x%(1<<i)+g%(1<<i))>=1<<i)<<(k-1-i) for i in range(1,k))
            assert rev(y,k)-rev(x,k)==rev(g,k)-3*c
            found.add(c)
        assert menu_size(k,g)==len(found)
        best=max(best,len(found))
    assert best==fib(k)
print('All carry languages through width 9 replayed, including the Fibonacci maximum.')
