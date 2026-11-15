def need(ok):
    if not ok:raise RuntimeError("matching decoder mismatch")
def rows(m):
    f=[0,1]
    for _ in range(m+1):f.append(f[-1]+f[-2])
    return f
def inverse(m,r,f):
    out=[]
    while m>0:
        if r>=f[m+1]:
            out.append(m)
            r-=f[m+1]
            m-=2
        else:m-=1
    need(r==0)
    return tuple(reversed(out))
for m in range(19):
    f=rows(m)
    a=set()
    for r in range(f[m+2]):
        p=inverse(m,r,f)
        need(p not in a and all(j-i>=2 for i,j in zip(p,p[1:])))
        need(r==sum(f[i+1] for i in p))
        a.add(p)
print("all finite Fibonacci matching sums through 18 switches")
