from math import factorial
# quotient order x D y S b
edges={(0,1),(0,3),(1,2),(1,3),(1,4),(2,0),(2,3),(3,4),(4,0),(4,2)}
def skeletons(q,i,j):
    tags=[0,2,4]+[1]*i+[3]*j
    n=len(tags);full=(1<<n)-1
    ends=[[0]*n for _ in range(1<<n)]
    for k in range(n):ends[1<<k][k]=1
    for mask in range(1,full+1):
        for last,h in enumerate(ends[mask]):
            if h:
                for u in range(n):
                    if not(mask>>u&1) and (tags[last],tags[u]) in edges:
                        ends[mask|1<<u][u]+=h
    H=[1]+list(map(sum,ends[1:]))
    v=[0]*(1<<n);v[0]=1
    for _ in range(q):
        out=[0]*(1<<n)
        for mask in range(1<<n):
            sub=mask
            while True:
                out[mask]+=H[sub]*v[mask^sub]
                if not sub:break
                sub=(sub-1)&mask
        v=out
    k=factorial(i)*factorial(j)
    assert v[full]%k==0
    return v[full]//k
for q in range(1,4):
    C=[[skeletons(q,i,j) for j in range(1,q+2)] for i in range(1,q+2)]
    print(q,C)

def Z(A):
    n=len(A);full=(1<<n)-1
    P=[[0]*n for _ in range(1<<n)]
    for i in range(n):P[1<<i][i]=1
    for m in range(1,full+1):
        for i,h in enumerate(P[m]):
            avail=A[i]&~m
            while h and avail:
                d=avail&-avail;P[m|d][d.bit_length()-1]+=h;avail-=d
    H=[1]+list(map(sum,P[1:]))
    Q=[0]*(1<<n)
    for m in range(1<<n):
        s=m
        while True:
            Q[m]+=H[s]*H[m^s]
            if not s:break
            s=(s-1)&m
    return H[full],Q[full],sum(H[s]*Q[full^s] for s in range(1<<n))

def build(D,S):
    tags=[0,2,4]+[1]*len(D)+[3]*len(S)
    n=len(tags); A=[0]*n
    for i,u in enumerate(tags):
        for j,v in enumerate(tags):
            if (u,v) in edges:A[i]|=1<<j
    for start,R in [(3,D),(3+len(D),S)]:
        for k,r in enumerate(R):A[start+k]|=r<<start
    assert all(((A[i]>>j)&1)+((A[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n))
    return A

assert [[skeletons(1,i,j) for j in (1,2)] for i in (1,2)]==[[13,6],[6,3]]
assert [[skeletons(2,i,j) for j in (1,2,3)] for i in (1,2,3)]==[[104,100,28],[100,92,26],[28,26,8]]
D=[(2,104,147,213,67,29,165,51),(110,184,210,212,1,156,178,17)]
assert [Z(A) for A in D]==[(345,5640,45225),(345,5640,44505)]
z=[Z(build(A,(0,))) for A in D]
assert [v[:2] for v in z]==[(34185,1352400),(34185,1332240)]
assert [h+6*q for h,q,r in z]==[8148585,8027625]
print("bridge skeletons and literal two-read collision PASS")
