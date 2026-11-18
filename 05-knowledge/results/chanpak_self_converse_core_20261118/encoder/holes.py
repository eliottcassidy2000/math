import itertools
from check import profile

def hcount(a):
    n=len(a);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            x=D[s][i]
            if not x:continue
            m=a[i]&~s
            while m:
                v=m&-m;m-=v;D[s|v][v.bit_length()-1]+=x
    return sum(D[-1])

def blue(n):
    q=[(i,j) for i in range(n) for j in range(i+2,n)]
    blocks=[{(i,j),(n-1-j,n-1-i)} for i,j in q if (i,j)<=(n-1-j,n-1-i)]
    base=[2<<i for i in range(n-1)]+[0]
    choices=[tuple((i,1<<j,j,1<<i) for i,j in orb) for orb in blocks]
    for bits in itertools.product([False,True],repeat=len(blocks)):
        a=base[:]
        for b,orb in zip(bits,choices):
            for i,u,j,v in orb:
                if b: a[i]|=u
                else:a[j]|=v
        yield a

def strong(a):
    n=len(a)
    # With increasing Hamilton frame any strong-component cut is a prefix.
    return all(any(a[v]&((1<<i)-1) for v in range(i,n)) for i in range(1,n))

if __name__=='__main__':
    universe={1}
    for n in range(3,9):
        S={}
        for a in blue(n):
            if strong(a):
                h=hcount(a)
                if h<75:S.setdefault(h,a)
        universe.update(S)
        print(n,sorted(S))
    D={i*j*j for i in universe for j in range(1,9,2) if j!=7 and i*j*j<75}
    print('SC omissions, by Busch and SCC factorization:', sorted(set(range(1,75,2))-D))
    a=(126,124,248,240,96,64,128,51)
    assert strong(a) and hcount(a)==63 and 63 not in D
