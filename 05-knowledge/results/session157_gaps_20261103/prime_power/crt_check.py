from math import gcd,lcm

def need(c, what):
    if not c: raise AssertionError(what)

def factor(n):
    d=2; out={}
    while d*d<=n:
        while n%d==0: out[d]=out.get(d,0)+1;n//=d
        d+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def hmod(n,m):
    if m==1:return 0
    return ((10*pow(4,n,9*m)+6*n-1)%(9*m))//9

def crt(a,A,b,B):
    if B==1:return a
    return a+A*((b-a)*pow(A,-1,B)%B)

def strong_n(m,r):
    # r must be odd if m even and a 3-adic unit if 3|m.
    fs=factor(m); L=1
    for p in fs:
        if p>=5:
            k=1
            while pow(4,k,p)!=1:k+=1
            L=lcm(L,k)
    powers=factor(L)
    for p,e in fs.items():powers[p]=max(powers.get(p,0), e-(p==2))
    R,Q=0,1
    for p,a in sorted(powers.items()):
        if not a: continue
        e=fs.get(p,0)
        precision=e-(p==2) if e else 0
        W=p**precision
        if precision:
            target=r%(p**e)
            candidates=[]
            for x in range(W):
                n=crt(R,Q,x,W)
                # Add a multiple of all prior coordinates and this one.
                n+=Q*W*(max(e,2)+1)
                if hmod(n,p**e)==target:candidates.append(x)
            need(bool(candidates),('no local solution',m,r,p,R,Q))
            x=candidates[0]
        else:x=0
        # Extend the selected lower p-adic digits arbitrarily (zero upper digits).
        Z=p**a
        R=crt(R,Q,x,Z);Q*=Z
    n=R+Q*(max(2,m.bit_length())+1)
    need(n>=1 and hmod(n,m)==r%m,('strong CRT',m,r,n,hmod(n,m)))
    return n

def tournament_count_mod(m,r):
    fm=factor(m); b=fm.get(3,0); q=3**b
    if not b: return hmod(strong_n(m,r),m)
    a=0
    while a<b and r%3**(a+1)==0:a+=1
    c=3**a
    u=0 if a==b else (r//c)%(3**(b-a))
    m0=m//q; rm0=(r*pow(c,-1,m0))%m0 if m0>1 else 0
    m3=3**(b-a)
    target=crt(rm0,m0,u,m3)
    modulus=m0*m3
    n=strong_n(modulus,target)
    return (c*hmod(n,m))%m

for m in range(1,251):
    for r in range(m):
        if m%2==0 and r%2==0:continue
        got=tournament_count_mod(m,r)
        need(got==r,('full result',m,r,got))
print('Constructive CRT controls: every eligible residue mod M=1..250')
for m in (385, 455, 1001, 3**5, 2**9*3**3, 7*11*13):
    for r in range(m):
        if m%2==0 and r%2==0:continue
        got=tournament_count_mod(m,r)
        need(got==r,('extra modulus',m,r,got))
print('Additional coupled and prime-power moduli all pass')
