"""Independent direct-matrix Smith audit. No producer imports or external packages."""
from itertools import product, combinations
from math import gcd
import json

def require(v,msg):
    if not v: raise RuntimeError(msg)

def smith(matrix):
    a=[row[:] for row in matrix]; n=len(a); k=0
    while k<n:
        nz=[(abs(a[i][j]),i,j) for i in range(k,n) for j in range(k,n) if a[i][j]]
        if not nz:break
        _,i,j=min(nz); a[k],a[i]=a[i],a[k]
        for row in a:row[k],row[j]=row[j],row[k]
        while True:
            changed=False
            for i in range(k+1,n):
                if a[i][k]:
                    q=a[i][k]//a[k][k]
                    a[i]=[v-q*w for v,w in zip(a[i],a[k])]
                    if a[i][k]:a[k],a[i]=a[i],a[k]
                    changed=True;break
            if changed:continue
            for j in range(k+1,n):
                if a[k][j]:
                    q=a[k][j]//a[k][k]
                    for i in range(n):a[i][j]-=q*a[i][k]
                    if a[k][j]:
                        for row in a:row[k],row[j]=row[j],row[k]
                    changed=True;break
            if changed:continue
            bad=next(((i,j) for i in range(k+1,n) for j in range(k+1,n) if a[i][j]%a[k][k]),None)
            if bad:
                i,_=bad;a[k]=[v+w for v,w in zip(a[k],a[i])];continue
            break
        if a[k][k]<0:a[k]=[-v for v in a[k]]
        k+=1
    ds=[a[i][i] for i in range(n)]
    require(all(a[i][j]==0 for i in range(n) for j in range(n) if i!=j),'Smith not diagonal')
    require(all(not ds[i] or ds[i+1]%ds[i]==0 for i in range(n-1)),'Smith divisibility')
    return ds

def det(a):
    a=[row[:] for row in a]; n=len(a)
    if not n:return 1
    prev=1; sign=1
    for k in range(n-1):
        if not a[k][k]:
            j=next((j for j in range(k+1,n) if a[j][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=a[i][j]*pivot-a[i][k]*a[k][j]
                require(num%prev==0,'nonexact Bareiss')
                a[i][j]=num//prev
        for i in range(k+1,n):a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]

def minor_factors(a):
    n=len(a); ds=[]; last=1
    for k in range(1,n+1):
        g=0
        for rr in combinations(range(n),k):
            for cc in combinations(range(n),k):
                g=gcd(g,abs(det([[a[i][j] for j in cc] for i in rr])))
        if g:
            require(last and g%last==0,'minor divisibility');ds.append(g//last)
        else:ds.append(0)
        last=g
    return ds

def matrix(n,u,plus):
    h=(n-1)//2;a=[[0]*h for _ in range(h)]
    for x in range(1,h+1):
        a[x-1][x-1]+=1
        y=(u*x)%n;s=1
        if y>h:y=n-y;s=-1
        a[x-1][y-1]+=s if plus else -s
    return a

def order(u,n):
    z=u%n;o=1
    while z!=1:z=z*u%n;o+=1
    return o

def prime_divisors(n):
    ans=[];p=2
    while p*p<=n:
        if n%p==0:
            ans.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:ans.append(n)
    return ans

def v2(n):
    v=0
    while n%2==0:v+=1;n//=2
    return v

def prediction(n,u,plus):
    ones=twos=zeros=cycles=0
    for d in range(3,n+1,2):
        if n%d:continue
        o=order(u,d);powers={pow(u,j,d) for j in range(o)}
        negative=d-1 in powers
        length=o//2 if negative else o
        count=sum(gcd(j,d)==1 for j in range(1,d))//(2*length)
        eps=-1 if negative else 1
        eta=eps*((-1)**length if plus else 1)
        ones+=count*(length-1);twos+=count*(eta==-1);zeros+=count*(eta==1);cycles+=count
    return [1]*ones+[2]*twos+[0]*zeros,cycles

cycle_cases=0
for n in range(1,8):
    for signs in product((-1,1),repeat=n):
        for plus in (False,True):
            a=[[0]*n for _ in range(n)]
            for i,s in enumerate(signs):a[i][i]+=1;a[i][(i+1)%n]+=s if plus else -s
            eps=1
            for s in signs:eps*=s
            eta=eps*((-1)**n if plus else 1)
            require(smith(a)==[1]*(n-1)+[1-eta],('signed cycle',n,signs,plus))
            cycle_cases+=1

matrices=minors=mod_vectors=sign_words=0
hostile={}
for n in range(3,66,2):
    for u in range(1,n):
        if gcd(u,n)!=1:continue
        profile=[v2(order(u,p)) for p in prime_divisors(n)]
        uniform=len(set(profile))==1
        for plus in (False,True):
            a=matrix(n,u,plus);ds=smith(a);pred,C=prediction(n,u,plus)
            require(ds==pred,('arithmetic Smith',n,u,plus,ds,pred))
            nonsingular=not any(d==0 for d in ds)
            expected=uniform and ((profile[0]!=1) if plus else (profile[0]>0))
            require(nonsingular==expected,('profile',n,u,plus,profile))
            require(abs(det(a))==(0 if not nonsingular else 2**C),('determinant',n,u,plus))
            matrices+=1
            if n<=11:
                require(minor_factors(a)==ds,('minor factors',n,u,plus));minors+=1
            if n<=9:
                for q in (2,3,4):
                    count=0
                    for x in product(range(q),repeat=len(a)):
                        if all(sum(v*w for v,w in zip(row,x))%q==0 for row in a):count+=1
                        mod_vectors+=1
                    predicted=q**ds.count(0)*gcd(q,2)**ds.count(2)
                    require(count==predicted,('kernel',n,u,plus,q))
            if n<=15:
                count=0
                for x in product((-1,1),repeat=len(a)):
                    if all(sum(v*w for v,w in zip(row,x))==0 for row in a):count+=1
                    sign_words+=1
                allowed=all(t==(1 if plus else 0) for t in profile)
                require(count==(2**C if allowed else 0),('full signs',n,u,plus))
            if n==15 and u==2:hostile['plus' if plus else 'minus']=ds
print(json.dumps({'status':'PASS','signed_cycle_matrices':cycle_cases,'odd_residue_matrices':matrices,'determinantal_minor_crosschecks':minors,'modular_vectors':mod_vectors,'literal_sign_words':sign_words,'hostile_N15_u2':hostile,'producer_imports':False,'range_odd_N':[3,65],'scope':'independent finite matrix and sign controls; uniform proofs in audit'},indent=2,sort_keys=True))
