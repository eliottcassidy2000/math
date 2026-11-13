from math import lcm,gcd

def req(ok):
    if not ok:raise RuntimeError('replay mismatch')

def factors(n):
    out=[]
    p=2
    while p*p<=n:
        if n%p==0:
            e=0
            while n%p==0:n//=p;e+=1
            out.append((p,e))
        p+=1
    if n>1:out.append((n,1))
    return out

def matmul(A,B,m):
    return [[sum(A[i][k]*B[k][j] for k in range(3))%m for j in range(3)] for i in range(3)]

def W(n,a,b,m):
    # sum_{k<n}(n-k) a^(n-k-1)b^k, valid without division
    A=[[a,a,1],[0,a,1],[0,0,b]]
    Z=[[1,0,0],[0,1,0],[0,0,1]]
    while n:
        if n&1:Z=matmul(Z,A,m)
        n//=2;A=matmul(A,A,m)
    return Z[0][2]

def Q(d,n,z,m):
    ans=pow(z,d*n+1,m)
    for e in (-1,1):
        w=W(n,pow(z,d,m),pow(z+e,d,m),m)
        for j in range(1,d):ans+=w*pow(z,d-j,m)*pow(z+e,j-1,m)
    return ans%m

def lift(n,L,p,e,f,target):
    for k in range(e):
        step=L*p**k; mod=p**(k+1)
        choices=[n+step*j for j in range(p) if f(n+step*j,mod)==target%mod]
        req(len(choices)==1);n=choices[0]
    return n,L*p**e

def clock(p,e,primes):
    d=3 if p==5 else 2
    L=lcm(p-1,3,2,*[a for q in primes if q<p for a in (q,q-1)])
    n,L=lift(0,L,p,e,lambda n,m:Q(d,n,2,m)*pow(pow(2,d*n,m),-1,m)%m,0)
    for q in primes:
        if q<=p:continue
        L=lcm(L,q-1)
        req(gcd(L,q)==1)
        n=next(n+L*j for j in range(q) if Q(d,n+L*j,1,q)!=0)
        L*=q
    if n==0:n=L
    req(Q(d,n,2,p**e)==0 and Q(d,n,1,3)==1)
    req(all(Q(d,n,1,q)!=0 for q in primes))
    return d,n

def source(k,a,odd_primes,H):
    O=1
    for p,e in odd_primes:O*=p**e
    N=(1<<k)*3**a*O
    d=0
    while d<a and H%3**(d+1)==0:d+=1
    # Divide the target by 3^d modulo N/3^d (H itself is divisible).
    T=H//3**d
    t=a-d
    pad_factor_exp=0
    if t and T%3==2:
        pad_factor_exp=lcm(2,*[(p-1)*p**(e-1) for p,e in odd_primes])
    M=N//3**d
    if pad_factor_exp:
        h0=(pow(2,pad_factor_exp,M)+1)%M
        T=T*pow(h0,-1,M)%M
    L=1<<max(0,k-1)
    if k:
        mod=1<<k
        n=((9*T-5)*pow(3,-1,mod)%mod)//2
    else:n=0
    if t:
        while n%3:n+=L
        # Use the nonsingular n=0 (mod 3) branch.
        L*=3
        for j in range(1,t):
            mod=3**(j+1)
            choices=[n+L*i for i in range(3) if Q(2,n+L*i,1,mod)==T%mod]
            req(len(choices)==1);n=choices[0];L*=3
    for p,e in odd_primes:
        L=lcm(L,p-1)
        n,L=lift(n,L,p,e,lambda n,m:Q(2,n,1,m),T)
    if n<k+1:n+=L*((k+1-n+L-1)//L)
    req(Q(2,n,1,M)==T%M)
    return d,pad_factor_exp,n

def solve(k,a,odds,target):
    N=(1<<k)*3**a
    for p,e in odds:N*=p**e
    if N==1:return
    V=[clock(p,e,[q for q,_ in odds]) for p,e in odds]
    v=1;qv=pow(2,k if V else max(k,1),N)
    for d,n in V:v=v*Q(d,n,1,N)%N;qv=qv*Q(d,n,2,N)%N
    req(gcd(v,N)==1 and qv%((1<<k)*__import__('math').prod(p**e for p,e in odds))==0)
    H=(-target*pow(v,-1,N))%N
    d,m,n=source(k,a,odds,H)
    u=Q(2,n,1,N)*pow(3,d,N)%N
    qu=Q(2,n,2,N)*pow(12,d,N)%N
    pu=Q(2,n,3,N)*pow(33,d,N)%N
    if m:
        # A cyclic wrap of TT_m: z^(m+2)+z[(z+1)^m-(z-1)^m].
        u=u*(1+pow(2,m,N))%N
        qu=qu*(pow(2,m+2,N)+2*(pow(3,m,N)-1))%N
        pu=pu*(pow(3,m+2,N)+3*(pow(4,m,N)-pow(2,m,N)))%N
    req(u==H)
    L=a+1
    if a>d+1:
        md=3**(a-d-1)
        # Compute reduced u and qu independently: residues mod N are multiples of 3^d.
        uu=(u//3**d)%md;qq=(qu//3**d)%md
        req(gcd(uu*qq,3)==1)
        while (2*pow(2,L,md)*qq-uu)%md:L+=1
    qu=qu*pow(2,L,N)%N;pu=pu*pow(3,L,N)%N
    for c in (2,6):
        A=(-3*u+6*qu+c*pu)%N
        req(A%3**a==0 and (A*qv-u*v-target)%N==0)

for k,a,odds,ts in [(4,4,[(5,2),(7,2)],(1,3,7,21,27,101)),
                    (0,0,[(7,3),(13,1),(43,1)],(0,7,341)),
                    (3,3,[],(1,3,5,7,9,27)),(0,4,[],(0,1,2,3,4,5,11)),(5,0,[],(1,7,21)),(4,0,[(5,2)],(1,3,7,11))]:
    for t in ts:solve(k,a,odds,t)
print('simultaneous physical residue decoder controls passed')

def literal(d,n,z):
    cut=(1,)+(0,)*(d-1); cut=cut*n; N=d*n+1
    def arc(i,j):
        if i==0:return cut[j-1]
        if j==0:return not cut[i-1]
        return i<j
    rows=[sum(1<<j for j in range(N) if i!=j and arc(i,j)) for i in range(N)]
    dp=[[0]*N for _ in range(1<<N)]
    h=[1]+[0]*((1<<N)-1)
    for i in range(N):dp[1<<i][i]=1
    for m in range(1,1<<N):
        h[m]=sum(dp[m])
        for i,c in enumerate(dp[m]):
            if c:
                b=rows[i]&~m
                while b:
                    u=b&-b;b-=u;dp[m|u][u.bit_length()-1]+=c
    q=[1]+[0]*((1<<N)-1)
    for _ in range(z):
        out=[0]*(1<<N)
        for m in range(1<<N):
            s=m
            while True:
                out[m]+=q[s]*h[m^s]
                if not s:break
                s=(s-1)&m
        q=out
    return q[-1]
for d in (2,3):
    for n in (1,2,3):
        for z in (1,2,3):req(Q(d,n,z,1000000007)==literal(d,n,z)%1000000007)
print('division-free pattern transfer independently replayed on literal tournaments')
