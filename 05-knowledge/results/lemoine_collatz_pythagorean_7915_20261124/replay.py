from collections import Counter
from math import gcd, isqrt

def step(x):
    z = 3*x-1
    return z//(z & -z)

def iterate(x, count):
    for _ in range(count):
        x = step(x)
    return x

def path(x, cap=5000):
    seen = {}
    while x not in seen:
        assert len(seen) < cap
        seen[x] = len(seen)
        x = step(x)
    return list(seen), list(seen)[seen[x]:]

def children(u, v):
    return (2*u-v,u), (2*u+v,u), (u+2*v,v)

def trip(m, n):
    return m*m-n*n, 2*m*n, m*m+n*n

f = Counter()
total = 0
for m in range(2, 361):
    for n in range(1, m):
        if (m-n)%2 != 1 or gcd(m,n) != 1:
            continue
        a,b,c = trip(m,n)
        k = 3*c-2*(a+b)
        assert (c-3*k)**2+4*(a-b)**2 == 8*k*k
        assert k == (m-2*n)**2+n*n
        if (m,n) != (2,1):
            p = tuple(sorted((n,abs(m-2*n)),reverse=True))
            assert (m,n) in children(*p)
            assert trip(*p)[2] == k < c
        f[k] += 1
        total += 1
def expected(k):
    if k == 1:
        return 1
    count = 0
    for p in range(2, isqrt(k)+1):
        if k%p == 0:
            if p%4 != 1: return 0
            count += 1
            while k%p == 0: k //= p
    if k > 1:
        if k%4 != 1: return 0
        count += 1
    return 3*2**(count-1)
# (m-2n)^2+n²=k implies m<=sqrt(5k), so m<=360 is ample.
for k in range(1, 10001):
    assert f[k] == expected(k)
print("primitive Euclid pairs:",total,"| intercept fibers independently checked: 1..10000")
print("root children:", [trip(*x) for x in children(2,1)])
for x in (1,5,17,7,9,15,13,29,61,63,65):
    route,cycle = path(x)
    print(x,route,"period",len(cycle),"basin",min(cycle))
assert [min(path(x)[1]) for x in (7,9,15,13,17,29,61,63,65)] == [5,5,1,5,17,1,17,5,1]
clocks = []
for rep in (1,5,17):
    cycle = path(rep)[1]
    js = [((3*x-1)&-(3*x-1)).bit_length()-1 for x in cycle]
    clocks.append((len(js),sum(js),len(js)+sum(js)))
    print("cycle clock",rep,js,clocks[-1])
assert clocks == [(1,1,2),(2,3,5),(7,11,18)]
assert tuple(map(sum,zip(*clocks))) == (10,15,25)
for e in range(1,251):
    z = 3**e+1
    assert iterate(2**e+1,e) == z//(z & -z)
for j in range(1,251):
    z = 3**j-1
    assert iterate(4**j-1,j) == z//(z & -z)
    assert iterate(4**j+1,2*j) == (3**(2*j)+1)//2
print("prefix laws checked: 1..250")
def isprime(x):
    return x>=2 and not any(x%i==0 for i in range(2,isqrt(x)+1))
for N in (7,9,15,21):
    pq = [(N-2*q,q) for q in range(2,N//2+1) if isprime(q) and isprime(N-2*q)]
    print("Lemoine",N,pq)
assert [x for x in range(5,10000,2) if step(x)==x-4] == [15]
