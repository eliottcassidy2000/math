from itertools import permutations

def need(x, info):
    if not x: raise AssertionError(info)

def rows(n):
    # r=0, y=1, x_i=i+1
    N=2*n+2; a=[0]*N
    for i in range(1,N):
        for j in range(i+1,N): a[i] |= 1<<j
    a[0] |= 1<<1
    for i in range(1,2*n+1):
        if i%2: a[0] |= 1<<(i+1)
        else: a[i+1] |= 1
    return a

def strong(a):
    N=len(a)
    for v in range(N):
        R=1<<v
        while True:
            S=R
            for i in range(N):
                if R>>i&1:R |= a[i]
            if R==S:break
        if R != (1<<N)-1:return False
    return True

def dp(a):
    N=len(a); D=[[0]*N for _ in range(1<<N)]
    for i in range(N):D[1<<i][i]=1
    for mask in range(1,1<<N):
        for i in range(N):
            v=D[mask][i]
            if v:
                for j in range(N):
                    if not(mask>>j&1) and a[i]>>j&1:
                        D[mask|1<<j][j]+=v
    return sum(D[-1])

def h(n):return (10*4**n+6*n-1)//9

for n in range(1,5):
    a=rows(n); need(strong(a),('strong',n)); need(dp(a)==h(n),('count',n,dp(a),h(n)))
print('independent subset DP: n=1..4:',[h(n) for n in range(1,5)])
for p in (5,7,11,13,17):
    L=1
    while pow(4,L,p)!=1:L+=1
    for e in (1,2,3):
        M=p**e
        values={h(L*t)%M for t in range(1,M+1)}
        need(len(values)==M,('prime-power residues',p,e,len(values)))
print('odd prime controls: all residues for p=5,7,11,13,17 and e<=3')
for e in range(1,13):
    M=2**e
    values={h(n)%M for n in range(e,e+M)}
    need(values==set(range(1,M,2)),('2-power residues',e))
print('2-power controls: all odd residues for e<=12')
