"""Replay Gaussian, Farey and physical tournament controls."""
from fractions import Fraction as F
from random import Random

ROOTS = ((8,37,1,6,15,29),(12,33,2,6,15,29))
UNIT71 = ((6,8,106,17,103,11,43),(6,8,122,17,99,11,43))


def ck(ok,*what):
    if not ok:
        raise RuntimeError(what)


def rooted_words(a):
    q = [0]*4
    full = (1 << len(a))-2
    def go(S,last,j):
        if S==full and j:
            q[j]+=1
        for v in range(len(a)):
            if last!=-1 and not a[last] >> v & 1:
                continue
            if v==0:
                if j<3:
                    go(S,0,j+1)
            elif not S>>v&1:
                go(S|1<<v,v,j)
    go(0,-1,0)
    return q[1:]


def strong(a):
    for i in range(len(a)):
        S=1<<i
        for _ in a:
            for j in range(len(a)):
                if S>>j&1: S|=a[j]
        if S!=(1<<len(a))-1:
            return False
    return True


def blow(a, parts):
    starts = [0]
    for p in parts:
        starts.append(starts[-1]+len(p))
    return [v << starts[i] | sum(((1 << len(p))-1) << starts[j]
            for j,p in enumerate(parts) if a[i] >> j & 1)
            for i,b in enumerate(parts) for v in b]


def count(a):
    size = 1 << len(a)
    dp = [[0]*len(a) for _ in range(size)]
    for i in range(len(a)):
        dp[1 << i][i] = 1
    for S,row in enumerate(dp):
        for i,v in enumerate(row):
            if not v: continue
            A = a[i] & (size-1-S)
            while A:
                b = A & -A
                A -= b
                dp[S|b][b.bit_length()-1] += v
    return [1]+[sum(row) for row in dp[1:]]


def pair(a):
    h = count(a)
    return h[-1], sum(x*y for x,y in zip(h,reversed(h)))


def tt(n):
    return [sum(1 << j for j in range(i+1,n)) for i in range(n)]


def sample(n,r):
    a = [0]*n
    for i in range(n):
        for j in range(i+1,n):
            x,y = (i,j) if r.randrange(2) else (j,i)
            a[x] |= 1 << y
    return a


def atan_bounds(a,m):
    t=sum(F((-1)**i,(2*i+1)*a**(2*i+1)) for i in range(m))
    u=t+F((-1)**m,(2*m+1)*a**(2*m+1))
    return min(t,u),max(t,u)


def main():
    a,b=atan_bounds(5,5)
    c,d=atan_bounds(239,5)
    lo,hi=16*a-4*d,16*b-4*c
    ck(F(223,71)<F(333,106)<lo<=hi<F(355,113)<F(22,7))
    for t,e in zip(ROOTS,(17,19)):
        ck(strong(t),t)
        ck(rooted_words(t)==[e,9,0],t)
    r=Random(239)
    src=[tt(j) for j in range(1,10)]+[sample(5,r) for _ in range(50)]
    for x in src:
        hh,pp=pair(x)
        for m,sgn in zip(ROOTS,(-1,1)):
            host=blow(m,[x]+[[0]]*5)
            ck(count(host)[-1]==9*pp+sgn*hh,x,sgn)
    u,v=UNIT71
    ck(strong(u) and strong(v) and sum((a^b).bit_count()
       for a,b in zip(u,v))==2)
    ck(pair(u)[0]==71 and pair(v)[0]==73)
    for b in range(1,26,2):
        for k in range(3,12):
           for s in [-1,1]:
              N=(b<<k)+s
              ck((3*N-s)//2==3*(b<<(k-1))+s and (3*N-s)%4==2)
              ck((3*N+s)//4==3*(b<<(k-2))+s and (3*N+s)%8==4)
    print("PASS: rational pi bounds, rooted words, 118 literal lifted hosts,",
          "the one-arc pair 71/73 and signed Collatz identities")


if __name__=='__main__':
    main()
