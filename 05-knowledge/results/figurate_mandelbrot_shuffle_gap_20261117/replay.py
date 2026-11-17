from itertools import combinations, permutations
from random import Random
from math import comb

def paths(T, zone):
    return [w for w in permutations(zone)
            if all(T[a]>>b&1 for a,b in zip(w,w[1:]))]

def shuffle_count(T,a,b):
    n,m=len(a),len(b)
    # two ending types and counts of ordered prefixes.
    D={(1,0,0):1,(0,1,1):1}
    for i in range(n+1):
        for j in range(m+1):
            for k in range(2):
                h=D.get((i,j,k),0)
                if not h:continue
                u=a[i-1] if k==0 else b[j-1]
                for v,vv in [(0,a[i] if i<n else None),(1,b[j] if j<m else None)]:
                    if vv is not None and T[u]>>vv&1:
                        key=(i+(v==0),j+(v==1),v)
                        D[key]=D.get(key,0)+h
    return sum(D.get((n,m,k),0) for k in range(2))

def defect(T,path,other):
    return any(T[b]>>c&1 and T[c]>>a&1
               for a,b in zip(path,path[1:]) for c in other)

def check(T,A,B):
    a,b=paths(T,A),paths(T,B)
    h=len(paths(T,A+B))
    dA=sum(defect(T,w,B) for w in a)
    dB=sum(defect(T,w,A) for w in b)
    m=len(a)*len(b)
    s=0
    for x in a:
        for y in b:
            k=shuffle_count(T,x,y)
            assert (k==1)==(not defect(T,x,B) and not defect(T,y,A)),(T,x,y,k)
            s+=k
    assert h>=s>=m+dA*len(b)+dB*len(a)-dA*dB,(T,h,m,dA,dB)
    assert (dA==dB==0)==(h==m),(T,h,m,dA,dB)
    return h-m

def tournament(n,x):
    rows=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if x>>k&1:rows[i]|=1<<j
        else:rows[j]|=1<<i
    return rows

def mixed_shift(n,d):
    shift,k,t=0,0,0
    while k+shift<=n:
        t+=comb(n-shift,k)
        shift+=d[k%len(d)]
        k+=1
    return t

def verify():
    cases=0
    for n in range(2,6):
        for packed in range(1<<comb(n,2)):
            T=tournament(n,packed)
            for mask in range(1,(1<<n)-1,2):
                check(T,tuple(i for i in range(n) if mask>>i&1),
                        tuple(i for i in range(n) if not mask>>i&1))
                cases+=1
    r=Random(911)
    for n in [6,7,8]:
        for k in range(70):
            T=tournament(n,r.randrange(1<<comb(n,2)))
            perm=r.sample(range(n),n); cut=r.randrange(1,n)
            check(T,tuple(perm[:cut]),tuple(perm[cut:]))
            cases+=1
    a=[mixed_shift(n,(1,2)) for n in range(70)]
    b=[mixed_shift(n,(2,1)) for n in range(70)]
    for n in range(5,70):
        assert a[n]==2*a[n-1]-a[n-2]+a[n-5]
        assert b[n]==2*b[n-1]-b[n-2]+b[n-5]
    print('partition/shuffle cases',cases)
    print('offsets (1,2)',a[:15])
    print('offsets (2,1)',b[:15])
    print('PASS')
if __name__=='__main__':
    verify()

def prof(rows,r):
    n=len(rows)
    Q=(1<<n)-1
    P=[[0]*n for _ in range(Q+1)]
    H=[1]+[0]*Q
    for s in range(1,Q+1):
        for v in range(n):
            if not s>>v&1:continue
            u=s^(1<<v)
            P[s][v]=sum(P[u][j] for j in range(n) if rows[v]>>j&1) if u else 1
        H[s]=sum(P[s])
    E=sum(P[s][r]*H[Q^s] for s in range(1,Q+1) if s.bit_count()%2==0)
    O=sum(P[s][r]*H[Q^s] for s in range(1,Q+1) if s.bit_count()%2==1)
    return H[Q],H[Q^(1<<r)],E,O

def witness():
    for p,T,r,z in ((5,(2,4,1),0,410),(7,(2,36,9,19,7,29),1,529),
                    (11,(8,5,33,22,7,27),3,102827),(13,(2,4,25,19,35,15),2,2195059),
                    (17,(2,4,1),0,6995913425)):
        H,g,E,O=prof(T,r)
        if O!=E+g:raise ValueError(T)
        for bit in reversed(bin(z)[2:]):
            if bit=='0':E,O=E+O,E+O
            else:H,E,O=H+2*E,E+O-g,E+O
        if H!=(2**p+1)**2+2:raise ValueError((p,H))
        print("period",3*p,"quotient count",H,"center count",H*(2**(p-1)-1))

if __name__=='__main__':
    witness()
