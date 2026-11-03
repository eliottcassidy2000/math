"""Exact source-word inversion in the inherited 65-presentation seam bank."""
from functools import lru_cache
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
PREVIOUS=HERE.parent/'golden_seam_recurrence_20261103'
SELECTED=('s14','s85','s248','s261','s262','s268','s330','s453')
GOLDEN=('a','b','q0','q1','q2','q4','q5')
ATLAS=json.loads((PREVIOUS/'atom_atlas.json').read_text(encoding='utf-8'))
BANK={k:tuple(map(tuple,z['matrix'])) for k,z in sorted(ATLAS.items())
      if k in GOLDEN or k.startswith('r') or k in SELECTED}
COST={k:len(ATLAS[k]['source_rows']) for k in BANK}
ROOT=(1,0,0)


def need(ok,message):
    if not ok:raise RuntimeError(message)


def mm(A,B):
    return tuple(tuple(sum(a*b for a,b in zip(row,col)) for col in zip(*B)) for row in A)


def step(v,P):
    return tuple(sum(v[i]*P[i][j] for i in range(3)) for j in range(3))


def state(word):
    v=ROOT
    for c in reversed(word):v=step(v,BANK[c])
    return v


def adjugate(P):
    cof=[]
    for i in range(3):
        row=[]
        for j in range(3):
            a=[[P[r][s] for s in range(3) if s!=j] for r in range(3) if r!=i]
            row.append((-1)**(i+j)*(a[0][0]*a[1][1]-a[0][1]*a[1][0]))
        cof.append(row)
    det=sum(P[0][j]*cof[0][j] for j in range(3))
    return det,tuple(zip(*cof))


INVERSES={c:adjugate(P) for c,P in BANK.items()}


def defect(v):
    h,x,y=v
    return x+y-h,y-x,h+x-y


def valid(v):
    if v==ROOT:return True
    h,x,y=v
    return (h>=1 and h%2==1 and 0<=x<=16*h and 0<=y<=16*h
            and min(defect(v))>=0)


def predecessors(v):
    """Legal arithmetic candidates; reachability is checked recursively."""
    if not valid(v) or v==ROOT:return
    rank=v[1]+v[2]
    for c,(det,adj) in INVERSES.items():
        need(det!=0,('invertible bank',c))
        numerator=step(v,adj)
        if any(a%det for a in numerator):continue
        u=tuple(a//det for a in numerator)
        if valid(u) and u[1]+u[2]<rank:
            need(step(u,BANK[c])==v,('inverse multiplication',c,v,u))
            yield c,u


@lru_cache(None)
def polynomial(v):
    """(source vertex order, number of literal source words), all orders."""
    if v==ROOT:return ((0,1),)
    if not valid(v):return ()
    result={}
    for c,u in predecessors(v):
        for n,count in polynomial(u):
            m=n+COST[c]
            result[m]=result.get(m,0)+count
    return tuple(sorted(result.items()))


@lru_cache(None)
def witness(v):
    if v==ROOT:return ()
    if not valid(v):return None
    for c,u in predecessors(v):
        tail=witness(u)
        if tail is not None:return (c,)+tail
    return None


def forward_states(bound):
    seen={ROOT};stack=[ROOT]
    while stack:
        v=stack.pop()
        for P in BANK.values():
            h,x,y=v
            nh=h*P[0][0]+x*P[1][0]+y*P[2][0]
            if nh>bound:continue
            w=(nh,h*P[0][1]+x*P[1][1]+y*P[2][1],h*P[0][2]+x*P[1][2]+y*P[2][2])
            need(valid(w) and w[1]+w[2]>v[1]+v[2],('forward guard',v,w))
            if w not in seen:seen.add(w);stack.append(w)
    return seen


def forward_polynomials(states):
    """Independent bottom-up word count; no adjugates or inverse calls."""
    counts={ROOT:{0:1}}
    for v in sorted(states,key=lambda z:(z[1]+z[2],z)):
        for c,P in BANK.items():
            w=step(v,P)
            if w not in states:continue
            row=counts.setdefault(w,{})
            for n,m in counts[v].items():
                d=n+COST[c];row[d]=row.get(d,0)+m
    return {v:tuple(sorted(row.items())) for v,row in counts.items()}


def ternary_left(v,offset):
    need(offset in (2,4) and v[0]>=3,'ternary scope')
    h,x,y=v;c=offset//2
    if x==h+c:return ('a',)
    if x+y==h+c:return ('a','a')
    return None


def extension_search(v,target):
    """Independent saturation of every left attachment below target H."""
    seen={v};stack=[v];answer=None
    while stack:
        z=stack.pop()
        if z[0]==target:answer=z
        for P in BANK.values():
            w=step(z,P)
            if w[0]<=target and w not in seen:seen.add(w);stack.append(w)
    return answer,seen


def literal(word):
    rows=()
    for c in word:
        atom=ATLAS[c]['source_rows'];n,d=len(rows),len(atom)
        rows=tuple(a|(((1<<d)-1)<<n) for a in rows)+tuple(a<<n for a in atom)
    n=len(rows);out=[(((1<<n)-1)^a)&~(1<<i) for i,a in enumerate(rows)]
    for i in range(n-1):out[i]|=1<<(i+1);out[i+1]&=~(1<<i)
    return tuple(out)


def hamilton(rows):
    n=len(rows)
    if not n:return 1
    dp=[{} for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for end,count in dp[mask].items():
            candidates=rows[end]&~mask
            while candidates:
                bit=candidates&-candidates;candidates-=bit;j=bit.bit_length()-1
                row=dp[mask|bit];row[j]=row.get(j,0)+count
    return sum(dp[-1].values())


def endpoint_stats(rows):
    """H, end-at-0, and covers whose distinguished path starts at0."""
    n=len(rows);dp=[[0]*n for _ in range(1<<n)];start=[0]*(1<<n)
    for j in range(n):dp[1<<j][j]=1
    fixed=[[0]*n for _ in range(1<<n)];fixed[1][0]=1
    for mask in range(1,1<<n):
        for end in range(n):
            todo=rows[end]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1
                dp[mask|bit][j]+=dp[mask][end]
                fixed[mask|bit][j]+=fixed[mask][end]
        start[mask]=sum(fixed[mask])
    counts=[sum(row) for row in dp];counts[0]=1;full=(1<<n)-1
    cover=sum(start[mask]*counts[full^mask] for mask in range(1,1<<n) if mask&1)
    return counts[full],dp[full][0],cover


def words(n):
    if n==0:yield ()
    for c,size in COST.items():
        if size<=n:
            for tail in words(n-size):yield (c,)+tail


def low_defect_words():
    """Complete finite source-word list with d=1 or2, proved in THEORY."""
    candidates={('a',)*n for n in range(1,6)}
    candidates|={('a',)*i+(c,)+('a',)*j for c in BANK if c!='a'
                for i in range(3) for j in range(3)}
    return {w:state(w) for w in candidates if defect(state(w))[0] in (1,2)}
