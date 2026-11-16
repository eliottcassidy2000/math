def need(ok):
    if not ok:
        raise AssertionError
def step(n):
    return (3*n+1)//2 if n&1 else n//2
a=27; o=0
for k in range(1,71):
    o+=a&1
    a=step(a)
    if k<59:need(a>=27 and 3**o >=2**k)
    if k==59:need(a==23 and o==37 and 3**o<2**k)
    if k<70:need(a!=1)
need(a==1)
for m,want in enumerate((27,291,12439,2173995791,1449330527,966220351),1):
    p=3**m
    e=next(e for e in range(2*3**(m-1)) if (41*pow(2,e,p)+1)%p==0)
    a=2**m*((41*2**e+1)//p)-1
    need(a==want)
    for _ in range(m):
        need(a&1)
        a=step(a)
    need(a==41*2**e)
print("PASS: 27 first descent and the first six growing exact-tail examples")
