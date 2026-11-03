import itertools, math, json

def check(x,msg):
    if not x: raise RuntimeError(msg)

def mul_cd(i,j,n):
    if n==0:return (1,0)
    h=1<<(n-1); a,b=i//h,j//h; u,v=i%h,j%h
    conj=lambda x: 1 if x==0 else -1
    if a==0 and b==0:s,w=mul_cd(u,v,n-1);return s,w
    if a==0 and b==1:s,w=mul_cd(v,u,n-1);return s,w+h
    if a==1 and b==0:s,w=mul_cd(u,v,n-1);return s*conj(v),w+h
    s,w=mul_cd(v,u,n-1);return -s*conj(v),w

def comp(a,b):return tuple(a[b[i]] for i in range(len(a)))
def power(a,k):
    z=tuple(range(len(a)))
    for _ in range(k):z=comp(a,z)
    return z

gens=[]
for a in range(1,8):
    p=[]
    for k in range(16):
        s,j=mul_cd(a,k%8,3);p.append(j+8*((k//8)^(s<0)))
    gens.append(tuple(p))
identity=tuple(range(16));group={identity};stack=[identity]
while stack:
    g=stack.pop()
    for a in gens:
        h=comp(a,g)
        if h not in group:group.add(h);stack.append(h)
check(len(group)==128,'octonion group order')
pairs_tested=0
for p,q in [(2,3),(3,4),(3,5),(4,5),(5,7),(7,8)]:
    valid={(a,b) for a in group for b in group if power(a,p)==power(b,q)}
    predicted={(power(g,q),power(g,p)) for g in group}
    check(valid==predicted,'nilpotent torus-knot representation theorem')
    pairs_tested+=len(group)**2
s3=list(itertools.permutations(range(3)))
valid=[(a,b) for a in s3 for b in s3 if power(a,2)==power(b,3)]
check(len(valid)==12 and sum(comp(a,b)!=comp(b,a) for a,b in valid)==6,'S3 hostile control')

# Literal labelled tournaments: fixed 0->1; Gray cycle of the other arcs.
GRAY=[(0,1),(1,1),(1,0),(0,0)]
def tournament(phases):
    n=3*len(phases);a=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):a[i][j]=True
    for b,t in enumerate(phases):
        x,y=GRAY[t]
        for (i,j),v in [((3*b,3*b+2),x),((3*b+1,3*b+2),y)]:a[i][j]=bool(v);a[j][i]=not v
    return a

def hamilton(a):
    n=len(a);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if dp[mask][i]:
                for j in range(n):
                    if not mask>>j&1 and a[i][j]:dp[mask|1<<j][j]+=dp[mask][i]
    return sum(dp[-1])
def H(phase):return 3**phase.count(0)
state_count=0;edge_checks=0
for k in range(1,4):
    for phase in itertools.product(range(4),repeat=k):
        a=tournament(phase);check(hamilton(a)==H(phase),'literal Hamilton path count');state_count+=1
        for i in range(k):
            b=list(phase);b[i]=(b[i]+1)%4;b=tournament(b)
            check(sum(a[u][v]!=b[u][v] for u in range(3*k) for v in range(u+1,3*k))==1,'literal flip');edge_checks+=1

def formula(a,b):
    periods=[4//math.gcd(x,4) for x in a];L=math.lcm(*periods);total=0;k=len(a)
    for mask in range(1<<k):
        residues=[];ok=True
        for i in range(k):
            if mask>>i&1:
                g=math.gcd(a[i],4)
                if b[i]%g:ok=False;break
                m=4//g;r=0 if m==1 else (-b[i]//g)*pow(a[i]//g,-1,m)%m
                residues.append((r,m))
        if not ok:continue
        if any((r-s)%math.gcd(m,n) for r,m in residues for s,n in residues):continue
        total+=(1<<mask.bit_count())*(L//math.lcm(*(m for r,m in residues)))
    return total
orbit_checks=0
for k in range(1,4):
    for a in itertools.product(range(4),repeat=k):
        L=math.lcm(*(4//math.gcd(x,4) for x in a))
        for b in itertools.product(range(4),repeat=k):
            actual=sum(H(tuple((b[i]+t*a[i])%4 for i in range(k))) for t in range(L))
            check(actual==formula(a,b),'CRT Hamilton orbit sum');orbit_checks+=1

# Finite cover restriction: slope translation orbits = lifted components.
lift_checks=0
for d in range(1,13):
    for p in range(-8,9):
        for q in range(-8,9):
            if math.gcd(p,q)!=1:continue
            unseen=set(itertools.product(range(d),repeat=2));components=[]
            while unseen:
                start=next(iter(unseen));v=start;orb=[]
                while v not in orb:
                    orb.append(v);unseen.remove(v);v=((v[0]+p)%d,(v[1]+q)%d)
                components.append(orb)
            check(len(components)==d and all(len(o)==d for o in components),'slope lift degree/components');lift_checks+=1
print(json.dumps({'status':'PASS','octonion_group_order':len(group),'knot_representation_pairs':pairs_tested,'S3_trefoil_representations':12,'S3_noncommuting_pairs':6,'literal_tournament_states':state_count,'single_arc_checks':edge_checks,'CRT_orbit_checks':orbit_checks,'primitive_slope_cover_checks':lift_checks},indent=2,sort_keys=True))
