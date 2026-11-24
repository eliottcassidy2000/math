from itertools import combinations
from functools import cache

@cache
def height(n, reverse):
    rows = [(1<<n)-(1<<(i+1)) for i in range(n)]
    for i,j in reverse:
        rows[i]^=1<<j
        rows[j]^=1<<i
    f = [[0]*n for _ in range(1<<n)]
    for i in range(n):
        f[1<<i][i]=1
    for s,row in enumerate(f):
        for i,v in enumerate(row):
            if v:
                a=rows[i]&~s
                while a:
                    b=a&-a
                    a^=b
                    f[s|b][b.bit_length()-1]+=v
    return sum(f[-1])

for n in range(3,12):
    chords = [(a,b) for a,b in combinations(range(n),2) if b-a>=2]
    short = [(i,i+2) for i in range(n-2)]
    for e in chords:
        a,b=e
        d=b-a
        assert height(n,(e,))==1+2**(d-1)
        out={}
        for f in short:
            if e==f:continue
            i,j=f
            h=height(n,tuple(sorted((e,f))))
            out[i]=h
            if d==2:
                expect=5 if abs(i-a)==1 else 9
            elif i==a or j==b:
                expect=5 if d==3 else 3+3*2**(d-3)
            elif j<=a or b<=i:
                expect=3*(1+2**(d-1))
            else:
                expect=3+3*2**(d-2)
            assert h==expect,(n,e,f,h,expect)
        if d>=3:
            assert {k for k,v in out.items() if v==min(out.values())}=={a,b-2}
print("All one/two reversed-chord formulas checked through order eleven.")

for a in range(1,5):
    for b in range(1,5):
        n=a+b+1
        switch=tuple((min(a,j),max(a,j)) for j in range(n) if j!=a)
        assert height(n,switch)==(2**a-1)*(2**b-1)+2
print("One-point real projective switch checked on all 16 widths 1..4.")
