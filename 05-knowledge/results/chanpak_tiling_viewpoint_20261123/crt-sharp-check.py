from math import gcd

TABLE = {
 13:(1025,65,2,156,[9,146,151,142,103]),
 41:(257,1025,6,1640,[1262,1633,832,1151,203]),
 257:(65,257,3,65792,[37175,57720,62171,2839,29857]),
}

def primes(n):
    while n>1:
        p=next(i for i in range(2,n+1) if n%i==0)
        yield p
        while n%p==0: n//=p

def rank(rows,p):
    rows=[list(a) for a in rows]
    r=0
    for j in range(len(rows[0])):
        i=next((i for i in range(r,len(rows)) if rows[i][j]%p), None)
        if i is None: continue
        rows[i],rows[r]=rows[r],rows[i]
        t=pow(rows[r][j],-1,p)
        rows[r]=[(x*t)%p for x in rows[r]]
        for i in range(r+1,len(rows)):
            t=rows[i][j]
            rows[i]=[(x-t*y)%p for x,y in zip(rows[i],rows[r])]
        r+=1
    return r

for q,(u,v,g,d,logs) in TABLE.items():
    assert d==q*(q-1) and v%q==0 and v%(q*q) and gcd(v,33)==1
    assert gcd(77*129,q)==1
    units=[5,17,37,751,u]
    local=[(9,2,6,[5,3,0,2,3 if u==1025 else 5 if u==257 else 1]),
           (121,2,110,[74,49,42,38,61 if u==1025 else 52 if u==257 else 65]),
           (q*q,g,d,logs)]
    column = [[int(x%4==3) for x in units],[int(x%8 in (3,5)) for x in units]]
    orders = [2,2]
    for m,g,d,xs in local:
        assert pow(g,d,m)==1 and all(pow(g,d//p,m)!=1 for p in primes(d))
        assert all(pow(g,a,m)==x%m for x,a in zip(units,xs))
        column.append(xs)
        orders.append(d)
    receipt={}
    for p in primes(2*6*110*d):
        sub=[col for col,o in zip(column,orders) if o%p==0]
        got=rank(list(zip(*sub)),p)
        assert got==len(sub)
        receipt[p]=got
    print(q,receipt)
