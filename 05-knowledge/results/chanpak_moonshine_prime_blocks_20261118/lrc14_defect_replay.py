from math import lcm

CLOCKS = [
    (3,2,0),(5,4,3),(11,10,1),(23,11,8),(67,66,22),
    (199,99,76),(463,231,147),(3697,1848,1769),(4621,4620,3467),
    (49,21,13),
]
L = lcm(*(d for q,d,r in CLOCKS))
def require(p, info):
    if not p:
        raise AssertionError(info)
for q,d,r in CLOCKS:
    require(pow(2,d,q)==1 and all(pow(2,i,q)!=1 for i in range(1,d)),(q,d))
    require([i for i in range(d) if (38*pow(2,i,q)+1)%q==0]==[r],q)
S=[j for j in range(L) if all(j%d!=r for q,d,r in CLOCKS)]
require((L,len(S),sum(j%3==1 for j in S))==(27720,4716,1392), (L,len(S)))

def f(x):
    h=1
    while x:
        q=x//2
        if not x&1: h+=q+(q&1)
        x=q
    return h

def f1(b,e):
    return b*e+f(b-1)*f(e-1)-(b&1)*(e&1)

EXACT = [(29,2,8283070098),(33,2,133767957234),(37,2,2058347761263),
         (45,71,17552393816392)]
for j,b,e in EXACT:
    require(f1(b,e)==(38<<j)+1,(j,b,e))
    require(j in S,j)

# Below j=27 the old independently corrected 2**32-1 cap applies.
# Outside the residue wheel, a fixed divisor descends. At a surviving
# exponent, strip a single 7 and test the resulting odd core directly.
require((38<<26)+1 < 2**32 <= (38<<27)+1,'previous cap')
unresolved=[]
checks=[]
for j in range(27,366):
    if j%L not in S or any(j==k for k,b,e in EXACT):
        continue
    N=(38<<j)+1
    core=N//7 if N%7==0 else N
    witness=next((a for a in (2,3,5,7,11,13,17,19,23)
                  if pow(a,core-1,core)!=1),None)
    if witness is None:
        unresolved.append(j)
    else:
        checks.append((j,witness))
require(unresolved==[365],unresolved)
print('exact Fermat descents',len(checks),'bases',sorted({a for j,a in checks}))

# Independently count the small 49 factor.
rows=[0,5,17,39,75,87,143,63]
n=len(rows)
require(all((rows[i]>>j&1)+(rows[j]>>i&1)==1
            for i in range(n) for j in range(i)),rows)
dp=[[0]*n for _ in range(1<<n)]
for i in range(n): dp[1<<i][i]=1
for s,r in enumerate(dp):
    for i,h in enumerate(r):
        t=rows[i]&~s
        while t and h:
            b=t&-t
            t^=b
            dp[s|b][b.bit_length()-1]+=h
require(sum(dp[-1])==49,rows)

# Basic finite checks for the analytic relay bound; the proof is uniform.
for z in range(1<<12):
    require(1<=f(z)<=z+1,z)
    for m in range(1,4):
        require(f((z+1)*(1<<m)-1)==f(z),(z,m))
        if z:
            require(f(z<<m)==f(z)+z*((1<<m)-1)+(z&1),(z,m))
for a,g in [(13190,37),(10542,175)]:
    require(8192<a<a+g<16384,(a,g))
    require(19*512<a<a+g<19*1024,(a,g))
print("PASS: 27720-clock, four literal F1 endpoints, H49 and relay constants")
