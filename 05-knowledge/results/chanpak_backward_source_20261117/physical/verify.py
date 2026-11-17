def app(a,q):
 n=len(a);return tuple(x|((1-(q>>i&1))<<n) for i,x in enumerate(a))+(q,)
def profile(a):
 n=len(a);M=1<<n;d=[[0]*n for _ in range(M)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,M):
    for i,v in enumerate(d[s]):
        if not v:continue
        t=a[i]&~s
        while t:
            b=t&-t;t-=b;d[s|b][b.bit_length()-1]+=v
 h=[1]+[sum(r) for r in d[1:]]
 return h[-1],sum(h[t]*h[M-1-t] for t in range(M))
def ck(x):
 if not x:raise ArithmeticError('check failed')
a=(174,124,120,240,65,80,1,118)
b=list(a);b[5]^=64;b[6]^=32;b=tuple(b)
for u in (43,75,612,392):a,b=app(a,u),app(b,u)
ck(sum((x^y).bit_count() for x,y in zip(a,b))==2)
ck(profile(a)==(271787,6325688) and profile(b)==(271797,6325680))
ck(12*271787+15*6325688==98146764 and 3*271787+4*6325688==26118113)
ck(98146764-26118113==72028651)
print('PASS independent literal n12 profiles, edge and reader coefficients')
