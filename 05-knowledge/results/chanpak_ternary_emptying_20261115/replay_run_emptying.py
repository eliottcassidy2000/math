from itertools import product,combinations,permutations
from collections import Counter

def check(ok,*args):
    if not ok:raise AssertionError(args)
def count(w):
    return 1+sum(max(2,1<<(j-i-1)) for i,x in enumerate(w)
            for j in range(i+1,len(w)) if x=='1' and w[j]=='0')
seen=Counter()
w63=[]
for l in range(2,14):
    for z in product('01',repeat=l-2):
        w='1'+''.join(z)+'0'
        N=count(w)
        S=0
        for b in w:S=2*S+int(b)
        k=w.count('10')
        # Python str.count matches all occurrences here: "10" cannot overlap itself.
        check((N-k-1-int(S%3==2))%3==0,w)
        check(N not in (7,21),w)
        if N==63:w63.append(w)
        seen[(k,N%3)]+=1
check(w63==['1001100','1100110'],w63)
for m in range(1,8):
    t=3**m
    for e in (0,1):
        R={(pow(2,3*j,t)-pow(2,j,t)+3)%t for j in range(e,2*3**(m-1)+e,2)}
        check(R==set(range(0,t,3)),m,e)
for r in range(1,12):
    x=1<<r
    check(count('1'+'0'*r+'1'*r+'0'*r)==x**3-x+3,r)

def lit(rows):
    n=len(rows);full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for j in range(n):dp[1<<j][j]=1
    for s in range(1,full):
        for j,a in enumerate(dp[s]):
            if a:
                nxt=rows[j]&~s
                while nxt:
                    u=nxt&-nxt;nxt-=u;dp[s|u][u.bit_length()-1]+=a
    return sum(dp[full])

def inflate(rows,v):
    n=len(rows)
    old=[i for i in range(n) if i!=v]+[v]*3
    out=[0]*(n+2)
    for i,j in combinations(range(n+2),2):
        if old[i]==old[j]:
            direct= (j-i)==1
        else: direct=bool((rows[old[i]]>>old[j])&1)
        a,b=(i,j) if direct else (j,i)
        out[a]|=1<<b
    return out
rng_words=((2,4,1),(6,12,8,1),(6,12,24,17,3))
for rows in rng_words:
    for v in range(len(rows)):
        A=[0]*4
        V=[j for j in range(len(rows)) if j!=v]
        for j in (1,2,3):
            for p in permutations(V):
                for gaps in combinations(range(len(V)+1),j):
                    out=[]
                    for i in range(len(V)+1):
                        if i in gaps:out.append(v)
                        if i<len(V):out.append(p[i])
                    if all(rows[u]>>w&1 for u,w in zip(out,out[1:])):
                        A[j]+=1
        check(lit(inflate(rows,v))==3*A[1]+6*(A[2]+A[3]))
rows=(4,5,152,3,11,31,63,123)
D=[r.bit_count() for r in rows]
aut=0
for p in permutations(range(8)):
    if any(D[i]!=D[p[i]] for i in range(8)):continue
    if all(((rows[i]>>j)&1)==((rows[p[i]]>>p[j])&1) for i,j in combinations(range(8),2)):aut+=1
check(aut==1 and lit(rows)==63)
print('canonical apex-run controls:',sum(seen.values()),'63 words:',w63)
print('ternary precision through 3^7; checked physical triangle inflation and automorphism control')
print({w:count(w) for w in ('100','10010','1000110','1100110','1001100')})
