"""Independent exact arithmetic controls; no floating point needed."""
from math import gcd

def check(cond):
    if not cond:
        raise ValueError("Exact finite control failed")

def norm(x,y):
    return x*x+x*y+y*y

for k in range(1,9):
    q=1<<k
    out={norm(x,y)%q for x in range(q) for y in range(q)}
    got=[]
    for S in (0,1):
        d={(x*y+y*(S-x-y)+x*(S-x-y))%q for x in range(q) for y in range(q)}
        check(d=={(S*S*pow(3,-1,q)-a)%q for a in out})
        got.append(d)
    check(got[0]|got[1]==set(range(q)))
    want={0,*[a for a in range(1,q) if ((a&-a).bit_length()-1)%2==0]}
    check(out==want)
    check(out|{(a+1)%q for a in out}==set(range(q)))
    odd_readout = {(1+2*a)%q for a in out}
    check(odd_readout|{(a+2)%q for a in odd_readout}==set(range(1,q,2)))
    check({(x+2*x*x)%(q//2) for x in range(q//2)}==set(range(q//2)))
    for r in range(k+1):
        h=1<<r
        for j in range(k+1):
            p=q>>j
            # Sum over subgroup characters exactly, by counting pairs
            # with a multiple of p as their difference.
            counts=[0]*p
            for x in range(h):
                counts[x%p]+=1
            full=p*sum(x*x for x in counts)
            check(full-h*h==max(0,p*h-h*h))
    for a in range(q):
        for c in (1,6,16,432):
            g=gcd(a*c,q)
            # Orthogonal rows split into blocks of exact size g.
            period=q//g
            check(sum((a*c*d)%q==0 for d in range(q))==g)
            check(all(a*c*d%q==0 for d in range(0,q,period)))
check({norm(x,y)%8 for x in range(8) for y in range(8)}=={0,1,3,4,5,7})
print("Dyadic Fourier conductor controls and exact norm images passed through q=256")
