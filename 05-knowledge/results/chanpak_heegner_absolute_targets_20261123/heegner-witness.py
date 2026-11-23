"""Verify literal cores for the nearest allowed odd Ramanujan target."""
from math import isqrt, prod

CORES = {
    3677: (920,269,913,420,874,519,943,818,32,266),
    4787: (384,805,201,867,815,325,915,570,132,293),
    91433: (2464,2577,2875,2435,2601,714,3487,1302,2098,2505,2879,160),
    163129: (2688,3885,2385,2005,131,1501,2195,774,593,2164,981,1464),
    47: (354,352,355,455,463,344,256,359,0),
    211: (470,72,90,193,202,223,128,6,254),
    229: (44,45,136,0,463,92,15,363,111),
    431: (46,168,194,100,239,4,163,297,127),
    2143: (372,345,106,193,460,282,288,103,140),
    13907: (734,1496,1682,1540,136,991,1436,1544,669,82,817),
    3: (2,4,1),
}

def count(R):
    n=len(R); N=1<<n
    assert all(0<=r<N and not(r>>i&1) for i,r in enumerate(R))
    assert all((R[i]>>j&1)+(R[j]>>i&1)==1 for i in range(n) for j in range(i))
    d=[[0]*n for _ in range(N)]
    for i in range(n):
        d[1<<i][i]=1
    for s in range(1,N):
        for i,x in enumerate(d[s]):
            if not x:continue
            z=R[i]&~s
            while z:
                b=z&-z;z-=b
                d[s|b][b.bit_length()-1]+=x
    return sum(d[-1])

def ordinary_apex(w):
    n=len(w)
    R=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]+[0]
    for i,c in enumerate(w):
        if c=='I':R[i]|=1<<n
        else:R[-1]|=1<<i
    return R

def apex_formula(w):
    return 1+sum((1<<(j-i-1))+(j==i+1)
        for i,a in enumerate(w) for j,b in enumerate(w) if i<j and a+b=='OI')

def rows_for_ordinal(factors):
    rows=[]
    for f in factors:
        n=len(rows)
        rows=[r|(((1<<len(CORES[f]))-1)<<n) for r in rows]+[x<<n for x in CORES[f]]
    return rows

if __name__=='__main__':
    for f,rows in CORES.items():
        assert all(f%j for j in range(2,isqrt(f)+1))
        assert count(rows)==f
    c=640320**3+744
    left=[3677,4787,91433,163129]
    right=[3,3,47,211,229,431,2143,13907]
    for target,fs in ((c-1,left),(c+3,right)):
        assert prod(fs)==target
        print(target, 'order',len(rows_for_ordinal(fs)),'factors',fs)
    for a in range(1,7):
        for b in range(1,10):
            for c0 in range(1,5):
                x,y,z=2**a,2**b,2**c0
                w='O'+'I'*a+'O'+'I'*b+'O'+'I'*c0
                assert apex_formula(w)==4*x*y*z+(2*y+1)*(z-x)+1-y
    for a,num in ((1,121),(5,32761)):
        R=ordinary_apex('O'+'I'*a+'OIII'+'O'+'I'*a)
        assert count(R)==num
        print("Nagell",num,"order",len(R))
