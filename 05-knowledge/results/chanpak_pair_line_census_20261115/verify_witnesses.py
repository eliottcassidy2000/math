from functools import lru_cache
from score import mark_w

def decode(n,z):
    rows=[0]*n; bit=0
    for i in range(n):
        for j in range(i+1,n):
            if z>>bit&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
            bit+=1
    return rows

def count(rows):
    @lru_cache(None)
    def endings(s,last):
        if s==1<<last: return 1
        r=s^(1<<last); w=r&~rows[last]
        ans=0
        while w:
            i=(w&-w).bit_length()-1; w&=w-1
            ans+=endings(r,i)
        return ans
    @lru_cache(None)
    def H(s):
        if not s: return 1
        return sum(endings(s,v) for v in range(s.bit_length()) if s>>v&1)
    return H

def close(a,m):
    x=m.bit_length(); n=len(a)
    v=x+n
    mask=sum((m>>i&1)<<(x-1-i) for i in range(x))
    rows=[]
    for i in range(x):
        rows.append((((1<<v)-1)^((1<<(i+1))-1)) | ((1<<v) if not mask>>i&1 else 0))
    rows += [(z<<x)|(1<<v) for z in a]
    rows += [mask]
    return rows

witnesses=[
 (32825,5,674,334),
 (34111,1,0,18168),
 (36833,1,0,20588),
]
for N,n,z,m in witnesses:
    A=decode(n,z)
    H=count(A)
    full=(1<<n)-1
    h=H(full); q=sum(H(s)*H(full^s) for s in range(1<<n))
    a=close(A,m)
    x=count(a)((1<<len(a))-1)
    w=mark_w(m)
    assert x==N==h+m*q-w*h, (N,x,n,z,m,h,q,w)
    print(N,n,z,m,h,q,w,len(a))

for m in range(1,2048):
    k=m//2
    w=mark_w(k)+(k+(k&1) if m&1 else 0)
    assert w==mark_w(m)
print('mask recurrence: 2047 exact controls')
