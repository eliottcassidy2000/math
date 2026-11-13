from itertools import combinations
from collections import Counter

def table(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    for j in range(n): dp[1<<j][j]=1
    for S in range(1,1<<n):
        for j,x in enumerate(dp[S]):
            if x:
                rest=rows[j]&~S
                while rest:
                    b=rest&-rest; rest-=b
                    dp[S|b][b.bit_length()-1]+=x
    return dp

def conv(t):
    n=len(t)
    return [sum(1<<j for j in range(n) if t[j]>>i&1) for i in range(n)]

def dbl(t,a,center=False,bridge=True,centerflip=False):
    n=len(t); N=2*n+center
    R=[0]*N; cross=(1<<n)-1
    for i in range(n):
        R[i]=t[i]|cross<<n
        R[n+i]=conv(t)[i]<<n
    if center:
        for i in range(n): R[i]|=1<<(2*n)
        R[2*n]=cross<<n
    if bridge: R[a]^=1<<(n+a); R[n+a]|=1<<a
    if centerflip:
        R[a]^=1<<(2*n); R[2*n]^=(1<<a)|(1<<(n+a)); R[n+a]|=1<<(2*n)
    return R

def folded(T,m,center):
    # enumerate half paths as an independent fixed-point check
    count=0
    def visit(S,v):
        nonlocal count
        if S==(1<<m)-1:
            w=2*m if center else (v+m)%(2*m)
            count+=T[v]>>w&1
            return
        for w in range(2*m):
            if not S>>(w%m)&1 and (v<0 or T[v]>>w&1):
                visit(S|(1<<(w%m)),w)
    visit(0,-1)
    return count

def decks(n):
    pairs=list(combinations(range(n),2))
    for code in range(1<<len(pairs)):
        r=[0]*n
        for k,(i,j) in enumerate(pairs):
            if code>>k&1: i,j=j,i
            r[i]|=1<<j
        yield r

if __name__=='__main__':
    stats=Counter()
    for n in range(1,5):
        for T in decks(n):
            D=table(T); C=table(conv(T)); N=(1<<n)-1
            h=sum(D[N])
            for a in range(n):
                e=D[N][a]
                L=sum(C[S][a]*(sum(D[N^S]) if N^S else 1)
                      for S in range(1,N+1) if S>>a&1)
                for c in (False,True):
                    t=dbl(T,a,c)
                    h1=sum(table(t)[-1]); f1=folded(t,n,c)
                    if h1!=(h*h+2*L*L if c else h*h+L*L-e*e): raise ValueError('bridge')
                    if f1!=(h if c else h+L-e): raise ValueError('fold')
                    stats[c]+=1
                m=sum(D[N^(1<<a)]) if n>1 else 1
                for bridge in (False,True):
                    t=dbl(T,a,True,bridge,True)
                    val=(h-e)**2+2*h*L+L*L+(2*L*L-2*m*(e+L) if bridge else 0)
                    if sum(table(t)[-1])!=val or folded(t,n,True)!=h-e+L:
                        raise ValueError((T,a,bridge,val,sum(table(t)[-1])))
                    stats['centerflip',bridge]+=1
                if L<h or L<e or ((L==e)!=(not T[a])): raise ValueError('defect')
    print(dict(stats))
    print('cycle input:',[(a,sum(table(dbl([2,4,1],a,c))[-1])) for a in range(1) for c in (False,True)])
