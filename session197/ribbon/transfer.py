"""Independent finite frontier transfer for span-w tournaments.
Indices 0..w-1 are the retained vertices (oldest first); -1 is forgotten.
Every component is (start,end). At most one start may be forgotten.
At each extension mask bit i=1 means new vertex points to retained vertex i.
"""
from collections import deque, defaultdict
from functools import lru_cache

def transitions(state,w,mask):
    pred=[None]+[i for i,(s,e) in enumerate(state) if e<0 or not(mask>>e&1)]
    succ=[None]+[i for i,(s,e) in enumerate(state) if s>=0 and mask>>s&1]
    ans=defaultdict(int)
    for a in pred:
        for b in succ:
            if a is not None and a==b: continue
            start=state[a][0] if a is not None else w
            end=state[b][1] if b is not None else w
            pairs=[p for i,p in enumerate(state) if i!=a and i!=b]+[(start,end)]
            shifted=[(s-1 if s>=0 else -1,e-1 if e>=0 else -1) for s,e in pairs]
            if sum(s<0 for s,e in shifted)>1: continue
            ans[tuple(sorted(shifted))]+=1
    return tuple(ans.items())

def automaton(w):
    states=[()]; idx={():0}; mats=[]
    for st in states:
        row=[]
        for mask in range(1<<w):
            ts=[]
            for nxt,c in transitions(st,w,mask):
                if nxt not in idx: idx[nxt]=len(states); states.append(nxt)
                ts.append((idx[nxt],c))
            row.append(ts)
        mats.append(row)
    return states,mats

def advance(vec,mat,mask,cap=None):
    out=[0]*len(vec)
    for i,x in enumerate(vec):
        if x:
            for j,c in mat[i][mask]:
                out[j]+=x*c
                if cap is not None and out[j]>cap: out[j]=cap
    return tuple(out)

def count(word,w):
    states,mat=automaton(w); v=(1,)+(0,)*(len(states)-1)
    for mask in word:v=advance(v,mat,mask)
    return sum(x for s,x in zip(states,v) if len(s)==1)

def graph(word,w):
    n=len(word); a=[[False]*n for _ in range(n)]
    for j,mask in enumerate(word):
        for i in range(j):
            backward=i>=j-w and (mask>>(i-(j-w)))&1
            a[i][j]=not backward;a[j][i]=bool(backward)
    return a

def held_karp(a):
    n=len(a)
    if not n:return 0
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for bits in range(1,1<<n):
        for i in range(n):
            if dp[bits][i]:
                for j in range(n):
                    if not(bits>>j&1) and a[i][j]:dp[bits|1<<j][j]+=dp[bits][i]
    return sum(dp[-1])

if __name__=='__main__':
    for w in range(1,5):
        states,mat=automaton(w)
        print(w,len(states),sum(len(t) for row in mat for t in row),flush=True)
