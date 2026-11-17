from math import gcd,isqrt
from itertools import combinations,product

AR=(158,164,56,242,98,1,167,52)
BR=(10,52,81,70,105,141,163,31)
A=((389,5600,1168,1376,102,88,217),(389,5600,1168,1376,102,88,225))
B=((521,7376,1662,1662,124,124,331),(521,7376,1662,1662,124,124,267))

def require(condition):
    if not condition:raise ArithmeticError('control failed')

def ends(a):
    n=len(a);N=1<<n;tab=[[0]*n for _ in range(N)]
    for i in range(n):tab[1<<i][i]=1
    for s in range(1,N):
        for i,v in enumerate(tab[s]):
            t=a[i]&~s
            while v and t:
                b=t&-t;t-=b;tab[s|b][b.bit_length()-1]+=v
    return tab

def p(a,tags):
    ends1=ends(a);n=len(a);N=1<<n
    sta=ends(tuple(sum(((a[j]>>i)&1)<<j for j in range(n)) for i in range(n)))
    H=[1]+[sum(row) for row in ends1[1:]]
    E,S=[[sum(x for j,x in enumerate(row) if tags>>j&1) for row in tab] for tab in (ends1,sta)]
    M=N-1
    P,W,Z,U=[sum(x[s]*y[M^s] for s in range(N)) for x,y in ((H,H),(S,H),(E,H),(E,S))]
    return H[M],P,W,Z,E[M],S[M],H[M]-U

for rows,states,marks in [(AR,A,(144,6)),(BR,B,(40,20))]:
    for i,r in enumerate(rows):
        z=1<<i
        for _ in rows:
            for j,x in enumerate(rows):
                if z>>j&1:z|=x
        require(z==(1<<len(rows))-1)
    require(tuple(p(rows,t) for t in marks)==states)

def score(word,k):
    J,W,E,K=1,0,0,1
    for h,p0,w,z,e,s,b in word:
        J,W,E,K=h*J,p0*W+J*w,J*e,h*K-J*(h-b)-E*s-W*z
    return W*2**k+K,J,W

def grid(j,M,k):
    L=(1<<(j-3))*M
    R=L if L%2 else L+M
    base,J,W=score([A[0]]*(2*L)+[B[0]]*(2*R)+[A[0]]*(L-1),k)
    unit=8*J//389; mod=8*L
    require(base%8==1 and gcd(unit//8,L)==1)
    z=(1-base)//8*pow(unit//8,-1,L)%L if L>1 else 0
    c=base+z*unit+L*unit-64*L*J//521
    a,b=L*unit,-64*L*J//521
    positions=[[(1,0),(-1,-1),(0,1)],[(-1,1),(0,0),(1,-1)],[(0,-1),(1,1),(-1,0)]]
    matrix=[[c+x*a+y*b for x,y in row] for row in positions]
    for row in positions:
        for x,y in row:
            val=score([A[1]]*((x+1)*L)+[A[0]]*((1-x)*L)+[B[1]]*((y+1)*L)+[B[0]]*(2*R-(y+1)*L)+[A[1]]*z+[A[0]]*(L-1-z),k)[0]
            require(val==c+x*a+y*b and val%mod==1)
    require(len({v for row in matrix for v in row})==9)
    require(all(sum(line)==3*c for line in matrix+list(zip(*matrix))+[[matrix[t][t] for t in range(3)],[matrix[t][2-t] for t in range(3)]]))
    return c,a,b

require(grid(3,1,3)==(2967778437612897,844724392,-5045647424))
for j,M,k in product(range(3,7),(1,3,5,7,15,35),range(3,8)):grid(j,M,k)
print('PASS: subset-DP tile profiles, fixed grid and 120 full modular grid replays')
