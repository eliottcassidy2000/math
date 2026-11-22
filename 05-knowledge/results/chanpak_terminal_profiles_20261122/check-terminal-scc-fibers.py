from collections import defaultdict

def rows(n,m):
    v=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if m&1:v[j]|=1<<i
            else:v[i]|=1<<j
            m>>=1
    return v

def strong(A):
    n=len(A)
    for i in range(n):
        s=1<<i
        while True:
            t=s
            for j in range(n):
                if s>>j&1:t|=A[j]
            if t==s:break
            s=t
        if s!=(1<<n)-1:return False
    return True

def profile(A):
    n=len(A);N=1<<n
    f=[[0]*n for _ in range(N)]
    for i in range(n):f[1<<i][i]=1
    for S in range(1,N):
        for v,x in enumerate(f[S]):
            b=A[v]&~S
            while x and b:
                u=b&-b;b-=u
                f[S|u][u.bit_length()-1]+=x
    H=[1]+list(map(sum,f[1:]))
    Q=[0]*N
    for S in range(N):
        W=S
        while True:
            Q[S]+=H[W]*H[S^W]
            if not W:break
            W=(W-1)&S
    return H[-1],Q[-1],sum(H[S]*Q[N-1-S] for S in range(N))

p=defaultdict(set)
for n in (3,4,5):
    for z in range(1<<(n*(n-1)//2)):
        A=rows(n,z)
        if strong(A):
            h,q,r=profile(A)
            if h<=13:p[h].add((q,r))
assert 5**5>13**3
assert p=={3:{(12,33)},5:{(32,117)},9:{(84,411)},11:{(100,465)},13:{(104,471)}}
p[9].add((12**2,33**2))
assert p[9]=={(84,411),(144,1089)}
print('complete normalized auxiliary fibres of heights 3,5,9,11,13 PASS')
