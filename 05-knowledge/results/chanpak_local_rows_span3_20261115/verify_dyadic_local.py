"""Exact controls: cancellation by the selected partner bit and primitive aliases."""
import random

P = [0]*16384
for z in range(2,len(P)):
    u, bit = divmod(z,2)
    P[z] = P[u]+(1-bit)*(u+(u%2))
assert all(v%2==0 for v in P)

def F(x,y):
    return 100*x*y+89*x+11*(P[x]+P[y]+1+x%2)

def G(x,y):
    return (F(x,y)-1)//2

def pairs(start,L,w):
    return [(i,i+w) for i in range(start,start+L)
            if (i-start)%(2*w)<w and i+w<start+L]

rng=random.Random(20261115)
for d in range(1,9):
    w=d+1
    for z in range(600):
        X0=[rng.randrange(len(P)),rng.randrange(len(P))]
        X1=X0[:]
        upper=[0,0]
        lower=[[],[]]
        for tape in range(2):
            start=rng.randrange(2)
            L=rng.randrange(14-start+1)
            for i,j in pairs(start,L,w):
                lower[tape].append(i)
                upper[tape]|=1<<j
                if rng.randrange(2): X1[tape]^=1<<i
        if X0==X1:
            continue
        tape=next(i for i in range(2) if X0[i]!=X1[i])
        bit=(X0[tape]^X1[tape]).bit_length()-1+w
        assert upper[tape] & (1<<bit)
        D0=G(*X1)-G(*X0)
        X0[tape]^=1<<bit
        X1[tape]^=1<<bit
        D1=G(*X1)-G(*X0)
        assert (D1-D0)%(2**d)==2**(d-1)

for j in range(16):
    k=j//2
    h=j-k
    pairs2=[(a*2**k,b*2**h-1) for a,b in ((3,15),(5,9))]
    assert all(F(x,y)==4500*2**j+77 for x,y in pairs2)

for a in range(1,1000,2):
    u=a//2
    g=P[a-1]//2
    q0=(a-1-P[a])//2
    assert g+q0==(3*u+u%2)//2
print("PASS: two-tape selected-bit checks, sixteen central primitive aliases, transpose shadow")
