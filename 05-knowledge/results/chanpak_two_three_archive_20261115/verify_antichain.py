def decode(n):
    r=[]
    a=0
    while n:
        if n&1:
            b=0;t=1
            while 3*t<=n:t*=3;b+=1
            r.append((a,b));n-=t
        n//=2;a+=1
    return r
def need(ok):
    if not ok:raise RuntimeError("decoder mismatch")
seen={}
for n in range(1,100000):
    r=decode(n)
    need(sum(2**a*3**b for a,b in r)==n)
    need(all(a<c and b>d for (a,b),(c,d) in zip(r,r[1:])))
    seen.setdefault(len(r),n)
for m,n in seen.items():need(n==3**m-2**m)
need(decode(19)==[(0,2),(1,1),(2,0)])
print("99,999 antichain decodes, exact first occurrences of greedy length")
