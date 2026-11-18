

def tt(n):return [(1<<n)-(1<<(i+1)) for i in range(n)]
def conv(a):
    return [sum((a[j]>>i&1)<<j for j in range(len(a))) for i in range(len(a))]
def ordinal(*parts):
    n=sum(map(len,parts));a=[];off=0
    for b in parts:
        e=off+len(b)
        a += [(r<<off) | ((1<<n)-(1<<e)) for r in b]
        off=e
    return a

def mirror(X,B):
    a=ordinal([0],X,B,conv(X),[0]); n=len(a)
    a[0]^=1<<(n-1);a[-1]^=1
    rho=[n-1]+list(range(len(X)+len(B)+1,n-1))+[len(X)+1+i for i in range(len(B)-1,-1,-1)]+list(range(1,len(X)+1))+[0]
    return a,rho

def profile(a):
    n=len(a);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            if not D[s][i]:continue
            m=a[i]&~s
            while m:
                v=m&-m;m-=v
                D[s|v][v.bit_length()-1]+=D[s][i]
    h=[1]+[sum(r) for r in D[1:]]
    return h[-1],sum(h[s]*h[((1<<n)-1)^s] for s in range(1<<n))

def orientations(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for w in range(1<<len(pairs)):
        a=[0]*n
        for k,(i,j) in enumerate(pairs):
            if w>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
        yield a

if __name__=='__main__':
    # 3-cycle relabeled so the anti-involution on its indices is reversal.
    c=[2,4,1]
    tests=0
    for n in range(1,5):
      for X in orientations(n):
        h,q=profile(X); reads=[]
        for B in ([],[0],c,tt(3)):
          b,p=profile(B);a,rho=mirror(X,B)
          assert all((a[i]>>j&1)==(a[rho[j]]>>rho[i]&1) for i in range(len(a)) for j in range(len(a)))
          g,_=profile(a)
          assert g==b*h*h+p*q*q
          reads.append(g);tests+=1
        g0,g1,g3,gt=reads
        assert (2*g0-g1,g1-g0)==(h*h,q*q)
        assert g3==3*h*h+12*q*q
        assert (2*g3-3*gt,3*gt-g3)==(3*h*h,12*q*q)
    print(tests,"literal wrapped counts and blue reflections through source order 4")
