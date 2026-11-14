#!/usr/bin/env python3
from functools import cache
from math import isqrt

@cache
def p(t):
    if not t:
        return 0
    u, b = divmod(t, 2)
    return p(u) + (0 if b else u + u % 2)

def rev(s, r):
    return int(f"{s:0{r}b}"[::-1], 2)

def decoded(y, h, S):
    if y <= 0:
        return [0] if y == 0 else []
    r = max(1, isqrt(y).bit_length())
    q = 1 << r
    out = []
    for s in range(q):
        a = q*S + h*(q-1-rev(s,r))
        d = 1-(s>>(r-1))
        m = y-S*s-h*p(s)
        if m < 0:
            continue
        lo = max(0, -(-(m-h*d)//(a+h)))
        hi = m//a
        for u in range(lo,hi+1):
            if a*u + h*(p(u)+d*(u%2)) == m:
                out.append(q*u+s)
    return sorted(out)

def count(rows):
    n = len(rows)
    @cache
    def f(k,x):
        if k == 1<<x:
            return 1
        return sum(f(k^(1<<x),v) for v in range(n) if (k & rows[x])>>v&1)
    return sum(f((1<<n)-1,x) for x in range(n))

def apex(rows, bits):
    r = list(rows) + [sum(bit << i for i,bit in enumerate(bits))]
    for j,bit in enumerate(bits):
        if not bit:
            r[j] |= 1 << len(rows)
    return r

def add_suffix(rows,k):
    n = len(rows)
    ext = ((1<<k)-1)<<n
    return [z|ext for z in rows]+[ext&~((1<<(i+1))-1) for i in range(n,n+k)]

checked=0
# Arbitrary strong/non-strong cores, and several independently chosen source masks.
for raw in ([],[0],[2,4,1],[8,1,3,6]):
    n = len(raw)
    rows = [((1<<(n+2))-2)] + [(z<<1)|(1<<(n+1)) for z in raw] + [0]
    h = count(rows)
    for cut in range(1<<n):
        mark = [1]+[(cut>>i)&1 for i in range(n)]+[0]
        n0 = count(apex(rows, mark))
        S = count(apex(add_suffix(rows,1), mark+[0]))-n0
        for t in range(1,65):
            bits = [1-(t>>i&1) for i in range(t.bit_length())]
            observed = count(apex(add_suffix(rows,len(bits)),mark+bits))
            expected = n0+S*t+h*p(t)
            if observed != expected:
                raise ValueError((rows,mark,t,observed,expected))
            checked += 1
print(f"literal Hamilton DP checks: {checked}")
for h,S in [(1,2),(5,64),(5,90),(31,341),(7,2)]:
    expect = {}
    for t in range(1001):
        y = S*t+h*p(t)
        expect.setdefault(y,[]).append(t)
    for y in range(2001):
        if decoded(y,h,S) != expect.get(y,[]):
            raise ValueError((y,h,S,decoded(y,h,S),expect.get(y,[])))
print("exact inverse checks: 10005")

for B in range(2, 81, 2):
    data = [(x,z) for x in range(B) for z in range(0,B,2)]
    def nextstates(x,z):
        return ((2*x%B,(z+x+x%2)%B),((2*x+1)%B,z))
    incoming = {s:[] for s in data}
    for s in data:
        for v in nextstates(*s):
            incoming[v].append(s)
    if any(len(w)!=2 for w in incoming.values()):
        raise ValueError((B,"not bistochastic"))
    for nxt in (lambda x,z: nextstates(x,z), lambda x,z: incoming[(x,z)]):
        todo, seen = [(0,0)], {(0,0)}
        while todo:
            for v in nxt(*todo.pop()):
                if v not in seen:
                    todo.append(v)
                    seen.add(v)
        if len(seen) != len(data):
            raise ValueError((B,"nonirreducible",len(seen),len(data)))
print("stationarity and irreducibility of digit chains: 40 even moduli")
