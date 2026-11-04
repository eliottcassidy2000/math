import itertools, json

def require(b, msg):
    if not b: raise RuntimeError(msg)
def graph(n, code):
    a=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        a[i][j]=bool(code>>k&1); a[j][i]=not a[i][j]
    return a
def hp(a, vs):
    return sum(all(a[x][y] for x,y in zip(p,p[1:])) for p in itertools.permutations(vs))
def profile(a):
    n=len(a); hs={mask:hp(a,[i for i in range(n) if mask>>i&1]) for mask in range(1,1<<n)}
    # Enumerate labelled assignments directly, not a cover convolution.
    q=[0]*(n+1)
    for k in range(1,n+1):
        for labels in itertools.product(range(k), repeat=n):
            masks=[sum(1<<i for i in range(n) if labels[i]==j) for j in range(k)]
            if not all(masks): continue
            z=1
            for m in masks:z*=hs[m]
            q[k]+=z
    return q
checks=0
for n in range(1,6):
    base=[0]*(n+1)
    for k in range(1,n+1):base[k]=sum(len(set(x))==k for x in itertools.product(range(k),repeat=n))
    for code in range(1<<(n*(n-1)//2)):
        a=graph(n,code); q=profile(a)
        require(all(q[k]>=base[k] for k in range(1,n+1)), ('bound',n,code))
        trans=not any(a[i][j] and a[j][k] and a[k][i] for i,j,k in itertools.permutations(range(n),3))
        require((q==base)==trans, ('equality',n,code))
        if not trans:require(q[1]>=3 and q[2]+2*q[1]>=3*2**(n-1),('nontrans',n,code))
        if n==4:require((q[1],q[2]) in {(1,14),(3,18),(5,22)},('fourprofile',code,q))
        if n==5:require(trans or 3*q[1]+q[2]>45,('five exclusion',code,q))
        checks+=1
print(json.dumps({'status':'PASS','labelled_profiles_orders_1_through_5':checks,'method':'literal permutations and labelled slot assignments','checks':['coefficientwise Stirling bound','transitivity equality iff','nontransitive F bound','order-four profiles','order-five low-spectrum exclusion']},indent=2))
