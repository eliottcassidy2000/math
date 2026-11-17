from growth import profile,add
from math import isqrt
from pathlib import Path
from collections import defaultdict
import re,sys

def read(f):
    g=a=None;B=defaultdict(dict)
    for l in Path(f).read_text().splitlines():
        if l.startswith('g='):g=int(l.split()[0][2:])
        q=re.match(r'(?:A (\d+)| d (-?\d+)) t=(\d+) h=(\d+) U=(\d+) D (.*)',l)
        if not q:continue
        a0,d,t,h,u,D=q.groups()
        if a0 is not None:a=int(a0)
        B[g,a][int(d or 0)]=(int(t),int(h),int(u),tuple(map(int,D.split())))
    return B

def witness(k,a,d,t,b):
    t0,h0,u,D=b
    n=len(D)
    A=tuple(x|((not (u>>i&1))<<n) for i,x in enumerate(D))+(u,)
    A=add(A,n,(1<<(t-t0))-1)
    j=len(A);far=[v for v in range(j) if v!=n]
    DD=tuple(sum((A[v]>>b&1)<<i for i,b in enumerate(far)) for v in far)
    U=sum((A[n]>>v&1)<<i for i,v in enumerate(far))
    g,_,E0,E,T0,T1,J=profile(DD,U)
    h=profile(A,0)[0]
    assert (h,g,2*E,g+E0-E)==(a*2**t-k+d,k,a*2**t,k)

if __name__=='__main__':
    B={}
    for f in sys.argv[1:]:B.update(read(f))
    P=defaultdict(set)
    for (g,a),b in B.items():
        for d in (-2,2):
            if d not in b:continue
            if (g,a)==(5,5):continue # all legal ready counts are divisible by 5 and >35
            for t in range(1,7):
                q=(a*2**t-g)%126
                if q%7==0 and q//7 in (1,5,7,11,13,17):P[q//7].add((g,a,t%6,d))
    assert set(P)=={5,7,11,13,17}
    print('phase support',sorted(P))
    for p,g,a,t,d in ((113,9,25,5,-2),(547,11,15,8,2),(173,5,19,6,2),(13,9,25,2,-2),(197,37,177,3,-2)):
        assert a*2**t-g==7*p and all(p%i for i in range(2,isqrt(p)+1))
        witness(g,a,0,t,B[g,a][0]);witness(g,a,d,t,B[g,a][d])
        print('prime residue / p / g / ready / companion:',p%18,p,g,7*p,7*p+d)
