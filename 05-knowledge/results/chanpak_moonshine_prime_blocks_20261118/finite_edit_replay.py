from itertools import product

def K(w):
    k=0
    while w:
        v=w//2
        if w&1:k+=v+(v&1)
        w=v
    return k

def out(p,w):
    h,e,t,d=p
    return h+t*w-d*K(w)-e*(w&1)

def decode(L,R,z):
    assert L[3]==R[3] and L[3]%2==1
    r=len(z)-1
    assert z[r]
    d=L[3]; D=R[2]-L[2]
    Z=sum(x*2**j for j,x in enumerate(z))
    A=2**(r+1)*D+d*sum(x*2**(r-j) for j,x in enumerate(z))
    assert A%2
    rhs=L[0]-R[0]+L[2]*Z-L[1]*z[0]
    sols=[]
    for lo in range(2**(r+2)):
        if any(x and ((lo>>j)&1)!=(x==-1) for j,x in enumerate(z)):continue
        B=D*lo+d*(K(lo+Z)-K(lo))+(L[1]-R[1])*(lo&1)
        num=rhs-B
        if num%(2*A):continue
        v=num//(2*A)
        if v<0:continue
        w=lo+2**(r+2)*v
        assert out(L,w+Z)==out(R,w)
        sols.append((w+Z,w,out(R,w)))
    return sols

profiles=[((1063,1,1847,1),(1063,1,1967,1)),((23,2,50,9),(25,2,50,9)),((17,5,40,3),(19,2,31,3))]
tested=0; examples=[]
for L,R in profiles:
 for r in range(4):
  for zs in product((-1,0,1),repeat=r+1):
   if not zs[-1]:continue
   z=list(zs); Z=sum(x*2**j for j,x in enumerate(z))
   dec=decode(L,R,z)
   brute=[]
   for w in range(4096):
    if any(x and ((w>>j)&1)!=(x==-1) for j,x in enumerate(z)):continue
    if out(L,w+Z)==out(R,w):brute.append((w+Z,w,out(R,w)))
   assert sorted(x for x in dec if x[1]<4096)==brute
   tested+=1
   if dec and L[0]==1063:examples.append((z,dec))
print({'mask_source_cases':tested,'exhaustive_addresses_per_case':4096,'physical_1063_examples':examples[:12], 'status':'PASS'})

def physical(w,r):
    A=[6,4,0]
    while w:
        bit=w&1; n=len(A)
        A=[(a|(1<<n))^(int(i==r)*bit<<n) for i,a in enumerate(A)]+[bit<<r]
        w//=2
    return A

def paths(A):
    n=len(A); dp=[[0]*n for _ in range(1<<n)]
    for j in range(n):dp[1<<j][j]=1
    for mask in range(1,1<<n):
        for j in range(n):
            if not dp[mask][j]:continue
            for k in range(n):
                if not (mask>>k&1) and (A[j]>>k&1):dp[mask|1<<k][k]+=dp[mask][j]
    return sum(dp[-1])
for z,expected in [([-1,0,1,-1,1],(22,11,37)),([-1,-1,-1,0,1],(24,15,47))]:
    assert expected in decode((1,0,2,1),(1,0,4,1),z)
    a,b,h=expected
    assert paths(physical(a,1))==paths(physical(b,0))==h
print('Literal multi-contact collisions at 37 and 47: PASS')
