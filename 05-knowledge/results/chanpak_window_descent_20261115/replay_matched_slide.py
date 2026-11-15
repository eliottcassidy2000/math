"""Reconstruct the physical core profiles, then replay the literal mark words."""
from fractions import Fraction
from replay_five_core import paths,starts,local,rev,require

R = (6,12,24,1,11)
u = 2
h = len(paths(R))
Q = h + sum(starts(R))
zs, ws = "101000010010", "010000000111"
d = len(zs)
L = 1 << d
z, w = int(zs,2), int(ws,2)
s, r = starts(R)[1], starts(rev(R))[1]
require((h,Q,s,r) == (13,104,18,18))
require((L-1)*s+h*z == Q*w)
require(Q*int(zs[::-1],2) == (L-1)*r+h*int(ws[::-1],2))
for i in range(len(R)):
    reachable = 1<<i
    for _ in R:
        reachable |= sum(1<<j for j in range(len(R)) if any(
            (reachable>>a&1) and (rows>>j&1) for a,rows in enumerate(R)))
    require(reachable == (1<<len(R))-1)

def profile(R, u):
    N,d = local(R,u)
    S = sum(v for i,v in enumerate(starts(R)) if u>>i&1)
    return N,S,d

def prepend_transitive(bits,R,u):
    n,m = len(bits),len(R)
    prefix = tuple((1<<(n+m))-(2<<i) for i in range(n))
    rows = prefix + tuple(x<<n for x in R)
    mark = (u<<n)|sum((c=="1")<<i for i,c in enumerate(bits))
    return rows,mark

def append(P, bits, h=h):
    N,S,d = P
    for bit in map(int,bits):
        if not bit:
            N += S+d
        S,d = 2*S+h*bit,h*bit
    return N,S,d

p,q = ws[-1],zs[0]
require(p==q=="1")
P0 = profile(*prepend_transitive(p,R,u))
require(P0 == (113,122,2))
require(append((0,0,0),zs) == (24076,33514,0))
require(append((0,1,0),zs) == (2938+24076,4096+33514,0))

for k in range(9):
    ref = None
    for j in range(k+1):
        x = int(p+ws*j,2)
        p1 = int((p+ws*j)[-1])
        s0 = Q*x+s
        d0 = local(R,u)[1]
        n0 = local(R,u)[0] + h*append((1,0,0),p+ws*j,h=1)[0] - h
        dt = local(rev(R), (1<<len(R))-1-u)[1]
        t = Q-h-r
        n0 += t*x + p1*dt
        N,S,D = append((n0,s0,d0),zs*(k-j)+q)
        carry = Fraction(9267808*4096**k+1944617,99225) + Fraction(9776*k,315)
        require(N == carry-20*j)
        yword = zs*(k-j)+q
        a,b = len(p+ws*j),len(yword)
        Tmark = (1<<a)*(Q*int(yword[::-1],2)+r)+h*int((p+ws*j)[::-1],2)
        both = (S,D,Tmark)
        if ref is None:
            ref = both
        require(ref == both)
        if k==1:
            require(N==(382625,382605)[j])
print("Strong 13-core; literal tile checks; profiles and integer counts match at all 45 test positions.")

def verify_case(R,u,zs,ws,expected):
    h = len(paths(R))
    Q = h + sum(starts(R))
    d = len(zs)
    s = sum(v for i,v in enumerate(starts(R)) if u>>i&1)
    r = sum(v for i,v in enumerate(starts(rev(R))) if u>>i&1)
    z,w = int(zs,2),int(ws,2)
    mask = (1<<len(R))-1
    require(((1<<d)-1)*s+h*z==Q*w)
    require(Q*int(zs[::-1],2)==((1<<d)-1)*r+h*int(ws[::-1],2))
    for i in range(len(R)):
        reach=1<<i
        for _ in R:
            for j,v in enumerate(R):
                if reach>>j&1: reach|=v
        require(reach==mask)
    Nc,ds = local(R,u)
    Nrc,de = local(rev(R),mask-u)
    require(Nc==Nrc)
    p,q = ws[-1],zs[0]
    def f(k,j):
        x,y = p+ws*j,zs*(k-j)+q
        a,b = len(x),len(y)
        sm = Q*int(x,2)+s
        Nm = Nc+h*(append((1,0,0),x,h=1)[0]-1)+(Q-h-r)*int(x,2)+int(x[-1])*de
        N,S,D = append((Nm,sm,ds),y,h=h)
        Tmark = (1<<a)*(Q*int(y[::-1],2)+r)+h*int(x[::-1],2)
        return N, (h,Q<<(a+b),S,Tmark,D,h*(1-int(x[0])))
    first = [f(1,j)[0] for j in range(2)]
    Delta = first[0]-first[1]
    require((h,abs(Delta))==expected)
    for k in range(5):
        coeffs = f(k,0)
        for j in range(k+1):
            v,profile = f(k,j)
            require((v,profile)==(coeffs[0]-j*Delta,coeffs[1]))
    print("additional literal slide h,Q,s,r",h,Q,s,r,"guard marks",p,q,
          "one-tile counts",first,"delta",Delta)
    return f,Delta

F13,D13=verify_case((6,12,24,1,11),2,"101001000010","010000001101",(13,10))
F29,D29=verify_case((70,8,2,5,15,31,62),49,"0001000011010010111100100010",
                    "1001010001000011000010110011",(29,6))
