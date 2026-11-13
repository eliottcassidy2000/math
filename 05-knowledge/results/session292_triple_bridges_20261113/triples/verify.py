#!/usr/bin/env python3
"""Independent DFS, squarefree pivot, divided-power controls; standard library."""
import itertools, math, random, json
from pathlib import Path

def ck(b,msg):
    if not b: raise RuntimeError(msg)
def add(*ps):
    q={}
    for p in ps:
        for s,a in p.items(): q[s]=q.get(s,0)+a
    return {s:a for s,a in q.items() if a}
def scale(p,a): return {s:a*b for s,b in p.items() if a*b}
def mul(p,q):
    out={}
    for s,a in p.items():
        for t,b in q.items():
            if not s&t: out[s|t]=out.get(s|t,0)+a*b
    return {s:a for s,a in out.items() if a}
def power(p,k):
    q={0:1}
    for _ in range(k): q=mul(q,p)
    return q
def mat(A):
    n=len(A); q=[[{} for _ in range(3*n)] for _ in range(3*n)]
    for i in range(n):
        q[n+i][i]={1<<i:1}; q[i][2*n+i]={0:1}
        for j in range(n):
            if A[i][j]: q[i][j]={1<<j:1}
    return q
def pivot(q,k,live):
    alive=[i for i in live if i!=k]
    for i in live:
        ck(all(s&(1<<k) for s in q[i][k]),'destination ideal')
        ck(not mul(q[i][k],q[k][k]),'loop annihilation')
    out=[[dict(p) for p in row] for row in q]
    for i in alive:
        for j in alive: out[i][j]=add(q[i][j],mul(q[i][k],q[k][j]))
    return out,alive
def dfs(A):
    n=len(A); out=[[{} for _ in A] for _ in A]
    def go(start,last,mask):
        out[start][last][mask]=out[start][last].get(mask,0)+1
        for nxt in range(n):
            if A[last][nxt] and not mask>>nxt&1: go(start,nxt,mask|1<<nxt)
    for i in range(n): go(i,i,1<<i)
    return out
rng=random.Random(292); tournaments=pivots=exchanges=triangle_action_checks=reciprocal_entries=0
for n in range(1,7):
    pairs=list(itertools.combinations(range(n),2)); size=1<<len(pairs)
    codes=range(size) if n<=4 else rng.sample(range(size),25)
    for code in codes:
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if code>>k&1: A[i][j]=1
            else: A[j][i]=1
        expected=dfs(A); orders=[list(range(n))]
        bar=[[{m:c*(-1)**m.bit_count() for m,c in cell.items()} for cell in row] for row in expected]
        brow=[add(*row) for row in bar]
        bcol=[add(*(bar[i][j] for i in range(n))) for j in range(n)]
        fbar=add({0:1},*brow); rr=add({0:1},scale(fbar,-1))
        ifbar=add(*(power(rr,k) for k in range(n+1)))
        ck(mul(fbar,ifbar)=={0:1},'reciprocal denominator')
        for i in range(n):
            for j in range(n):
                rhs=add(scale(bar[i][j],-1),mul(mul(brow[i],bcol[j]),ifbar))
                ck(expected[j][i]==rhs,'marked reciprocal identity')
                reciprocal_entries+=1
        identity=tuple(range(n))
        def compose(p,q): return tuple(p[q[i]] for i in range(n))
        def ppow(p,k):
            ans=identity
            for _ in range(k): ans=compose(ans,p)
            return ans
        pm=[]
        for perm in itertools.permutations(range(n)):
            if all(A[perm[i]][perm[j]]==A[i][j] for i in range(n) for j in range(n)) or all(A[perm[i]][perm[j]]==A[j][i] for i in range(n) for j in range(n)): pm.append(perm)
        for u in pm:
            if ppow(u,2)!=identity: continue
            for v in pm:
                if ppow(v,3)!=identity: continue
                for r in (5,7):
                    if ppow(compose(u,v),r)==identity: ck(u==v==identity,'triangle action')
                    triangle_action_checks+=1
        for _ in range(2):
            z=list(range(n)); rng.shuffle(z); orders.append(z)
        for order in orders:
            q=mat(A);live=list(range(3*n))
            for at,k in enumerate(order):
                if at+1<n:
                    u,v=k,order[at+1]
                    qu,lu=pivot(q,u,live); quv,luv=pivot(qu,v,lu)
                    qv,lv=pivot(q,v,live); qvu,lvu=pivot(qv,u,lv)
                    ck(luv==lvu,'live sets')
                    for i in luv:
                        for j in luv:
                            five=add(q[i][j],mul(q[i][u],q[u][j]),mul(q[i][v],q[v][j]),mul(mul(q[i][u],q[u][v]),q[v][j]),mul(mul(q[i][v],q[v][u]),q[u][j]))
                            ck(quv[i][j]==qvu[i][j]==five,'five term')
                    exchanges+=1
                q,live=pivot(q,k,live);pivots+=1
            ck([[q[n+i][2*n+j] for j in range(n)] for i in range(n)]==expected,'DFS-pivot')
        tournaments+=1
factorial_tests=0
for n in range(1,8):
    for _ in range(12):
        r={s:rng.randrange(-3,4) for s in range(1,1<<n) if rng.random()<.08}
        for k in range(n+2): ck(all(a%math.factorial(k)==0 for a in power(r,k).values()),'factorial divisibility')
        for a in range(4):
            for b in range(4):
                ra={s:c//math.factorial(a) for s,c in power(r,a).items()}
                rb={s:c//math.factorial(b) for s,c in power(r,b).items()}
                rab={s:c//math.factorial(a+b) for s,c in power(r,a+b).items()}
                ck(mul(ra,rb)==scale(rab,math.comb(a+b,a)),'divided power')
        factorial_tests+=1

def vp_fact(k,p):
    v=0
    while k: k//=p;v+=k
    return v
# Symmetric squarefree basis e_j(t_1,...,t_h); e_i*e_j=binom(i+j,i)e_(i+j).
def smul(a,b,h):
    return [sum(math.comb(k,i)*a[i]*b[k-i] for i in range(k+1)) for k in range(h+1)]
cutoffs=[]
for p in (2,3,5,7):
    for s in range(1,5):
        K=next(k for k in range(1,100) if vp_fact(k,p)>=s)
        h=K;D=[(-1)**i for i in range(h+1)]; r=[0]+[-D[i] for i in range(1,h+1)]
        q=[1]+[0]*h;inv=[0]*(h+1);prev=None
        for j in range(K+1):
            if j<K: inv=[a+b for a,b in zip(inv,q)]
            if j==K-1: prev=q[:]
            if j==K: ck(all(a%(p**s)==0 for a in q),'cutoff')
            q=smul(q,r,h)
        ck(prev[K-1]==math.factorial(K-1),'sharp coefficient')
        ck(prev[K-1]%(p**s)!=0,'strict sharpness')
        ck([a%(p**s) for a in smul(D,inv,h)]==[1]+[0]*h,'truncated inverse')
        cutoffs.append({'p':p,'s':s,'K':K,'triangles_for_sharpness':K-1})
# Rank-one single-use seam identity over the same occupancy ring.
def mmul(a,b):
    n=len(a)
    return [[add(*(mul(a[i][k],b[k][j]) for k in range(n))) for j in range(n)] for i in range(n)]
def resolvent(a):
    n=len(a); acc=[[{0:1} if i==j else {} for j in range(n)] for i in range(n)]
    term=[[dict(x) for x in row] for row in acc]
    for _ in range(n):
        term=mmul(term,a)
        acc=[[add(acc[i][j],term[i][j]) for j in range(n)] for i in range(n)]
    ck(not any(any(row) for row in mmul(term,a)),'resolvent nilpotence')
    return acc
seam_checks=0
for _ in range(30):
    n=4
    q=[[{1<<j:rng.randrange(1,4)} if rng.randrange(2) else {} for j in range(n)] for i in range(n)]
    a,b=rng.sample(range(n),2);delta={(1<<n)|(1<<b):2}
    qp=[[dict(x) for x in row] for row in q];qp[a][b]=add(qp[a][b],delta)
    r,rp=resolvent(q),resolvent(qp)
    for i in range(n):
        for j in range(n): ck(add(rp[i][j],scale(r[i][j],-1))==mul(delta,mul(r[i][a],r[b][j])),'single-use seam')
    seam_checks+=1
receipt={'status':'PASS' ,'tournaments':tournaments,'pivot_steps':pivots,'two_pivot_exchanges':exchanges,'random_divided_power_examples':factorial_tests,'triangle_action_checks':triangle_action_checks,'formal_seam_checks':seam_checks,'reciprocal_entries':reciprocal_entries,'cutoffs':cutoffs,'notes':['Exact DFS and polynomial path atlas agree.','Sharpness concerns geometric power, not number of cycles or complexity.']}
Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
