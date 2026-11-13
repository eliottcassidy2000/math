from itertools import permutations
from math import gcd

def req(t, msg):
    if not t: raise RuntimeError(msg)

def rows(n, code):
    a=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1: a[i]|=1<<j
            else: a[j]|=1<<i
            k+=1
    return a

def direct(a):
    n=len(a)
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for m in range(1,1<<n):
        for i,x in enumerate(d[m]):
            if not x: continue
            s=a[i]&~m
            while s:
                b=s&-s; s-=b; d[m|b][b.bit_length()-1]+=x
    return sum(d[-1])

def strong(a):
    m=1
    while True:
        z=m
        for i,b in enumerate(a):
            if m>>i&1:z|=b
        if m==z:break
        m=z
    if m!=(1<<len(a))-1:return False
    ar=[sum(1<<j for j,b in enumerate(a) if b>>i&1) for i in range(len(a))]
    m=1
    while True:
        z=m
        for i,b in enumerate(ar):
            if m>>i&1:z|=b
        if m==z:return m==(1<<len(a))-1
        m=z

def coeff(a,v,k):
    b=[i for i in range(len(a)) if i!=v]
    total=0
    for w in permutations(b+[v]*k):
        if all(a[x]>>y&1 for x,y in zip(w,w[1:])): total+=1
    from math import factorial
    return total//factorial(k)

def cycles(a,v):
    def go(x,m):
        s=int(m.bit_count()>=3 and a[x]>>v&1)
        q=a[x]&~m
        while q:
            b=q&-q;q-=b;s+=go(b.bit_length()-1,m|b)
        return s
    return go(v,1<<v)

def inflate(a,v,c):
    labels=[i for i in range(len(a)) if i!=v]+[v]*len(c)
    d=[]; q=len(a)-1
    for i,x in enumerate(labels):
        d.append(sum(1<<j for j,y in enumerate(labels) if
              (c[i-q]>>(j-q)&1 if x==y==v else a[x]>>y&1)))
    return d

for n in range(3,7):
    vals=set()
    for x in range(1<<(n*(n-1)//2)):
        a=rows(n,x)
        if strong(a): vals.add(direct(a))
    req(min(vals)=={3:3,4:5,5:9,6:15}[n] and 7 not in vals and 21 not in vals, "small spectrum")
    print(n,sorted(vals))
for n in range(2,5):
    for x in range(1<<(n*(n-1)//2)):
        a=rows(n,x)
        for v in range(n):
            c2,c3=(coeff(a,v,k) for k in (2,3))
            h=direct(a)
            r=inflate(a,v,[2,4,1]); s=inflate(a,v,[6,4,0])
            req(direct(r)==3*h+6*(c2+c3)==direct(s)+2*h,"substitution")
            if strong(a):
                req(c2>=cycles(a,v)>=1,"cycle injection")
src=[2,12,145,149,99,143,175,19]
req([direct(inflate(src,7,[(1<<m)-(1<<(i+1)) for i in range(m)])) for m in (1,2,3,4)]==[179,539,1619,4859],"ternary clock")
req(direct(inflate(src,7,[2,4,1]))==1977,"different triangle outcome")
print("substitution, cycle bound and literal clock pass")
for value,a in {
    35:[4,117,64,39,77,85,9], 49:[4,109,24,17,99,77,13],
    63:[94,88,194,20,4,95,24,123], 343:[226,68,49,231,139,146,52,70]
}.items():
    req(all(not(a[i]>>i&1) for i in range(len(a))) and all(
        (a[i]>>j&1)+(a[j]>>i&1)==1 for i in range(len(a)) for j in range(i+1,len(a))),
        "seed tournament")
    req(strong(a) and direct(a)==value,"exceptional seed")
print("literal strong seeds 35,49,63,343 pass")
