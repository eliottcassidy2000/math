from itertools import combinations

def need(t):
    if not t:raise AssertionError("guard replay")
def table(rows):
    n=len(rows);M=1<<n
    D=[[0]*n for _ in range(M)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,M):
        for i,v in enumerate(D[s]):
            if not v:continue
            c=rows[i]&~s
            while c:
                b=c&-c;c-=b;D[s|b][b.bit_length()-1]+=v
    return [1]+[sum(r) for r in D[1:]]
def guard(rows,open=()):
    n=len(rows);a=n;b=n+1
    return [r | ((1<<b) if i not in open else 0) for i,r in enumerate(rows)]+[(1<<n)-1,(1<<a)+sum(1<<i for i in open)]
def slots(h,k):
    D=[int(i==0) for i in range(len(h))]
    for _ in range(k):
        A=[]
        for s in range(len(h)):
            t=s;v=D[0]*h[s]
            while t:
                v+=D[t]*h[s^t]
                t=(t-1)&s
            A.append(v)
        D=A
    return D[-1]
values=set()
nlab=0
for n in range(1,6):
    edges=list(combinations(range(n),2))
    for s in range(1<<len(edges)):
        R=[0]*n
        for j,(a,b) in enumerate(edges):
            if s>>j&1:a,b=b,a
            R[a]|=1<<b
        B,G=table(R),table(guard(R))
        need(G[-1]==B[-1]+slots(B,2))
        if G[-1]<=63:values.add(G[-1])
        if s%13==0:
            for t in range(1,5):
                need(slots(G,t)==t*t*slots(B,t)+t*(slots(B,t+1)-slots(B,t-1)))
        nlab+=1
need(values=={3,5,9,15,17,27,33,37,51})
R=[((1<<6)-1)^((1<<(i+1))-1) for i in range(6)]
singles=[table(guard(R,(i,)))[-1] for i in range(6)]
need(singles==[95,79,69,61,51,33])
need([v for v in combinations(range(6),2) if table(guard(R,v))[-1]==63]==[(2,3)])
for i in range(1,25,2):
    # By linearity these monomials control arbitrary odd integral polynomials.
    f=[pow(j,i,9) for j in range(9)]
    for k in range(1,21):
        f=[(j*j*f[j]+j*(f[(j+1)%9]-f[j-1]))%9 for j in range(9)]
        need(f[1]%3==0)
        if k>=3:need((f[1]//3-(-1)**(k-3)*(i+(i==1)))%3==0)
print('guarded labelled hosts',nlab,'through 63:',sorted(values))
print('single contacts at n=6:',singles,'double witness (2,3):',table(guard(R,(2,3)))[-1])
print('240 mod-9 monomial iteration controls')
