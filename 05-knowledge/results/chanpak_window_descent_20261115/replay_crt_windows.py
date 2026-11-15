"""Arithmetic target certificates, and literal Hamilton checks on small cuts."""
from itertools import product
from math import gcd, prod, lcm
from functools import cache

def check(ok):
    if not ok:
        raise AssertionError("CRT-window control failed")

def decoder(t, hs, ms):
    p=prod(hs)
    ws=[p//h for h in hs]
    us=[t*pow(w,-1,h)%h for w,h in zip(ws,hs)]
    s=sum(w*u for w,u in zip(ws,us))
    carry, rem=divmod(t-s,p)
    check(rem==0)
    if not 0<=carry<=sum(ms)-len(hs):
        return None
    for i,(h,m) in enumerate(zip(hs,ms)):
        v=min(carry,m-1)
        us[i]+=h*v
        carry-=v
    check(carry==0 and sum(w*u for w,u in zip(ws,us))==t)
    return us

for hs,ms in [((3,5),(2,2)),((3,5),(1,2)),((3,5,7),(1,2,2)),
              ((5,7,9),(2,2,1)),((3,5,7),(1,1,1)),((3,7,11),(3,1,1))]:
    p=prod(hs)
    ws=[p//h for h in hs]
    k=len(hs)
    F=(k-1)*p-sum(ws)
    K=sum(ms)-k
    counts={sum(w*d for w,d in zip(ws,ds))
            for ds in product(*(range(h*m) for h,m in zip(hs,ms)))}
    for t in range(max(counts)+2):
        d=decoder(t,hs,ms)
        check((t in counts)==(d is not None))
        if d is not None: check(all(0<=x<h*m for x,h,m in zip(d,hs,ms)))
    check(F not in counts and (K+1)*p not in counts)
    check(set(range(F+1,(K+1)*p))<=counts)
for hs in [(3,3,5),(3,5,9),(15,21,35),(5,)*7]:
    check(gcd(*(prod(hs)//h for h in hs))==prod(hs)//lcm(*hs))
print("Exact target carries, both boundaries and repeated-catalogue gcd passed.")

def tt(n):
    return tuple((1<<n)-(2<<i) for i in range(n))

def join(a,b):
    return tuple(x|(((1<<len(b))-1)<<len(a)) for x in a)+tuple(x<<len(a) for x in b)

def apex(rows, mark):
    n=len(rows)
    return tuple(x | (int(not(mark>>i&1))<<n) for i,x in enumerate(rows))+(mark,)

@cache
def counts(rows):
    n=len(rows); f=(1<<n)-1
    A=[[0]*n for _ in range(f+1)]
    h=[1]+[0]*f
    for s in range(1,f+1):
        for i in range(n):
            if not s>>i&1: continue
            t=s^(1<<i)
            A[s][i]=sum(A[t][j] for j in range(n) if t&rows[i]&(1<<j)) if t else 1
        h[s]=sum(A[s])
    W=[sum(A[s][i]*h[s^f] for s in range(1,f+1)) for i in range(n)]
    return h[f],sum(h[s]*h[s^f] for s in range(f+1)),W

def conv(rows):
    return tuple(sum(1<<j for j,x in enumerate(rows) if x>>i&1)
                 for i in range(len(rows)))

def M(rows,x):
    return sum(((x>>(len(rows)-i-1))&1)<<i for i in range(len(rows)))

def rev(x,w):
    return M(range(w),x)

def rc(x,w):
    return ((1<<w)-1)^rev(x,w)

def padded_controls(core, pad, decode3=False):
    colors=range(1<<len(core))
    total=0
    table={}
    widest=[]
    for n in range(pad+1):
        groups={}
        for a in range(n+1):
            b=n-a
            R=join(join(tt(a),core),tt(b))
            h,Q,W=counts(R)
            E=counts(conv(R))[2]
            for x,y,U in product(range(1<<a),range(1<<b),colors):
                mark=M(range(a),x)|(U<<a)|(M(range(b),y)<<(a+len(core)))
                s=sum(w for i,w in enumerate(W) if mark>>i&1)
                t=sum(w for i,w in enumerate(E) if not(mark>>i&1))
                j=U.bit_count()
                hc=counts(apex(R,mark))[0]
                key=(h,Q,s,t)
                if decode3:
                    z=(x<<2)|j
                    z=(z<<b)|y
                    check((h,Q,s,t)==(3,12<<n,3*z,3*((4*rc(y,b)+3-j)*(1<<a)+rc(x,a))))
                    check(t//3-rc(z,n+2)==(rev(j,2)-j)*(1<<a))
                    if j in [0,3]:
                        check(hc==3*counts(apex(tt(n+2),M(range(n+2),z)))[0])
                groups.setdefault(key,set()).add(hc)
                total+=1
        widest.append(max(map(len,groups.values())))
        table.update(groups)
    return total,widest

n,w=padded_controls((2,4,1),3,True)
check(max(w)==1)
print("Count-three literal triangle-padding checks:",n)
n,w=padded_controls((6,12,8,1),3)
print("Count-five control: checks",n,"max (h,Q,S,T) widths by total padding 0..3",w)

def signed_end_profile(rows):
    n=len(rows); f=(1<<n)-1
    A=[[0]*n for _ in range(f+1)]
    for s in range(1,f+1):
        for i in range(n):
            if not s>>i&1: continue
            t=s^(1<<i)
            A[s][i]=sum(A[t][j] for j in range(n) if t&rows[i]&(1<<j)) if t else 1
    return A[f]

def signature(rows,mark):
    h,Q,W=counts(rows)
    flip=conv(rows)
    s=sum(w for i,w in enumerate(W) if mark>>i&1)
    t=sum(w for i,w in enumerate(counts(flip)[2]) if not(mark>>i&1))
    d=sum(w for i,w in enumerate(signed_end_profile(flip)) if mark>>i&1)
    f=sum(w for i,w in enumerate(signed_end_profile(rows)) if not(mark>>i&1))
    return h,Q,(s-d)//2,(s+d)//2,(t-f)//2,(t+f)//2

pairs=[((118,124,120,113,96,64,0),50,69),((126,124,120,48,96,64,8),38,67)]
for U,mark,hc in pairs:
    R=join(join((0,),U),(0,))
    m=(mark<<1)|(1<<(len(U)+1))
    check(signature(R,m)==(5,1024,102,107,149,154))
    check(counts(apex(U,mark))[0]==counts(apex(R,m))[0]==hc)
print("Literal count-five unit palette: matched response, c=31 or 32.")

# Intentionally nonconsecutive residues as integer sets, but consecutive lifts.
hs=(5,7)
palettes=([
    [u+5*v for v in range(u%2,3+(u==3))] for u in range(5)
],[
    [u+7*v for v in range(u%3-1,4-(u==2))] for u in range(7)
])
p=prod(hs)
ws=[p//h for h in hs]
A=sum(w*max(map(min,palette)) for w,palette in zip(ws,palettes))
B=sum(w*min(map(max,palette)) for w,palette in zip(ws,palettes))
reachable={sum(w*d for w,d in zip(ws,ds))
           for ds in product(*(sum(x,[]) for x in palettes))}
check(set(range(A-p+1,B+p))<=reachable)
check(A-p not in reachable and B+p not in reachable)
for z in range(-3,10):
    check((set(range(z*p,(z+1)*p))<=reachable)==(A//p<=z<=B//p))
print("Separable worst-lift formula passed on unequal carry ranges.")
