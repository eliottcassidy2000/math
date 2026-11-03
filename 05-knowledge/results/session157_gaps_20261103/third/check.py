
def need(x, info):
    if not x:raise AssertionError(info)

def word(n):return (1,)+(0,0,1)*n+(0,)

def tournament(n):
    w=word(n); N=len(w)+1; A=[0]*N
    for i in range(1,N):
        for j in range(i+1,N):A[i]|=1<<j
    for i,b in enumerate(w,1):
        if b:A[0]|=1<<i
        else:A[i]|=1
    return A

def strong(A):
    N=len(A)
    for v in range(N):
        R=1<<v
        while True:
            old=R
            for u in range(N):
                if R>>u&1:R|=A[u]
            if R==old:break
        if R!=(1<<N)-1:return False
    return True

def dp(A):
    N=len(A);D=[[0]*N for _ in range(1<<N)]
    for v in range(N):D[1<<v][v]=1
    for mask in range(1,1<<N):
        for i in range(N):
            c=D[mask][i]
            if c:
                for j in range(N):
                    if not(mask>>j&1) and A[i]>>j&1:D[mask|1<<j][j]+=c
    return sum(D[-1])

def F(n):
    x=80*8**n+28*n+67
    need(x%49==0,('integer',n))
    return x//49

def modF(n,m):return ((80*pow(8,n,49*m)+28*n+67)%(49*m))//49
for n in range(1,5):
    A=tournament(n);need(strong(A),('strong',n))
    need(dp(A)==F(n),('DP',n,dp(A),F(n)))
print('Independent graph DP n=1..4:',[F(n) for n in range(1,5)])
S={modF(4*t,75) for t in range(1,76)}
need(S==set(range(75)),('75 residues',sorted(S)))
print('F(4t), 1<=t<=75: all 75 distinct classes')
missing=[9,12,24,27,54,57,69,72]
for a in missing:
    t=next(t for t in range(1,76) if modF(4*t,75)==a)
    print('residue',a,'t',t,'n',4*t)
