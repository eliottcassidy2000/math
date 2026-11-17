from math import isqrt
from replay_mersenne import A as MINUS, need, strong, table, ordinal
PLUS=(8,1,3,6)

def tt(n):
    full=(1<<n)-1
    return tuple(full^((1<<(i+1))-1) for i in range(n))

def inflate(scaffold, block, terminal):
    labels=[(0,i) for i in range(len(block))]+[(i,0) for i in range(1,len(scaffold))]
    def arc(i,j):
        x,u=labels[i]; y,v=labels[j]
        return (block[u]>>v&1) if x==y==0 else (scaffold[x]>>y&1)
    return [sum(arc(i,j)<<j for j in range(len(labels))) for i in range(len(labels))],len(block)+terminal-1

def join(a,v,r):
    n=len(a)
    rows=list(a)+[((1<<n)-1)| (x<<n) for x in tt(r)]
    rows[v]|=1<<n
    rows[n]^=1<<v
    return rows

def profile(a,v):
    d,F=table(a)
    b=sum(row[v]*F[(len(F)-1)^s] for s,row in enumerate(d) if s>>v&1)
    need(d[-1][v]>0)
    return F[-1],b

def reads(block):
    _,f=table(block)
    full=len(f)-1
    return f[-1],sum(x*f[full^s] for s,x in enumerate(f))

def prime(x):
    return x>1 and all(x%d for d in range(2,isqrt(x)+1))

def verify(block,x,y,j=1,both=True):
    h,P=reads(block);b=(8**j-1)//7
    need(x==(1+4*b)*P+b*h and y==(16+48*b)*P-b*h)
    for setup,v,r,q in [(PLUS,1,3*j+1,x),(MINUS,5,3*j,y)]:
        a,m=inflate(setup,block,v)
        want=(3*P-h,4*P+h) if setup==PLUS else (64*P+h,96*P-2*h)
        need(profile(a,m)==want)
        if both:
            child=join(a,m,r)
            need(table(child)[1][-1]==7*q and strong(child))
    print(j,h,P,x,y)

for block,x,y in [(tt(1),11,127),(tt(7),641,8191),
                  (ordinal((24,1,3,6,14),tt(2)),2011,25589)]:
    need(prime(x) and prime(y))
    verify(block,x,y)

# Second repunit verifies changing right reader, not just the fixed four/three pair.
verify(tt(1),83,887,2)

# Cutoff for the negative tower at q=65537 (the singleton must be checked separately).
Q=65537
need(prime(Q))
checks = []
b=1
while 32+95*b<=Q:
    a=16+48*b
    need(2*a-b!=Q)
    M=4*a
    h=(-Q*pow(b,-1,M))%M
    hi=Q//(4*a-b) # every non-singleton exception with h>1 has order at least three
    # h=1 would require the residue to be one; each of these residues exceeds one.
    need(h>1)
    if h<=hi:
        need(h==255 and 7*a-b==447 and (7*a-b)*h>Q) # 5!<255 forces n>=6
    checks.append((b,h,hi))
    b=8*b+1
need(checks==[(1,255,257),(9,1479,36),(73,1031,4),(585,96519,0)])
print(checks)

for block,x,y in [
    (tt(2),21,255),
    ((24,1,3,6,14),511,6389),
    ((16,93,1,53,68,87,13),3011,38113),
]:
    verify(block,x,y)
    if y!=255:
        need(prime(y))
need({38113%8,4091%8,6389%8,8191%8}=={1,3,5,7})
need((1<<32)+1==641*6700417)
