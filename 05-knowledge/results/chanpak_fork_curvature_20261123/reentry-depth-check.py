"""Exact subset and finite-ring checks for the uniform reentry note."""
from math import factorial
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

spec = spec_from_file_location("restricted", Path(__file__).with_name("matching-restricted.py"))
lib = module_from_spec(spec)
spec.loader.exec_module(lib)

def profile(A):
    n=len(A)
    f=[[0]*n for _ in range(1<<n)]
    for i, row in enumerate(A):
        assert not row>>i&1
        assert all((row>>j&1)+(A[j]>>i&1)==1 for j in range(i))
        reach=1<<i
        for _ in A:
            for u,row in enumerate(A):
                if reach>>u&1:
                    reach|=row
        assert reach==(1<<n)-1
        f[1<<i][i]=1
    for S,row in enumerate(f):
        for i,c in enumerate(row):
            z=A[i]&~S
            while c and z:
                b=z&-z
                z-=b
                f[S|b][b.bit_length()-1]+=c
    h=[1]+[sum(r) for r in f[1:]]
    v=[1]+[0]*(len(h)-1)
    out=[]
    for _ in A:
        y=[0]*len(h)
        for S in range(1,len(h)):
            a=S
            while a:
                y[S]+=h[a]*v[S^a]
                a=(a-1)&S
        out.append(y[-1])
        v=y
    return out

A,B=lib.MA,lib.MB
PA,PB=profile(A),profile(B)
d=[b-a for a,b in zip(PA,PB)]
assert d==[0,0,240,480,240,0,0,0,0]
assert lib.direct(A)==lib.direct(B)==(2231,36088)
c=21*2231+840*36088
assert (c-2*2231,c,c+2*2231)==(30356309,30360771,30365233)
assert (lib.direct(lib.wrap(A)),lib.direct(lib.wrap(B)))==((38319,743964),(38319,744444))
assert [lib.direct(lib.wrap(lib.wrap(M)))[0] for M in (A,B)]==[782283,782763]
assert [m for m in (5,6,7,8,9) if all(x*factorial(i)%m==0 for i,x in enumerate(d))]==[5,6,8]
assert [lib.direct(lib.transitive(n)) for n in (3,4)]==[(1,8),(1,16)]
assert [lib.direct(lib.wrap(lib.transitive(n)))[1]%8 for n in (3,4)]==[4,0]

def mm(a,b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))%5
                 for i in range(2) for j in range(2))
for sign,expected in [(1,(6,24)),(-1,(3,12))]:
    a=(1,1,3,(4+2*sign)%5); I=(1,0,0,1)
    z=I; first_scalar=None
    for i in range(1,25):
        z=mm(z,a)
        if z[0]==z[3] and z[1:3]==(0,0) and not first_scalar:
            first_scalar=i
        if z==I:
            assert (first_scalar,i)==expected
            break
    else:
        raise AssertionError
print("Reentry twins, two successive sandwiches, and modular periods checked.")
