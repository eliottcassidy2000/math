#!/usr/bin/env python3
"""Independent literal-permutation and three-run-assignment two-apex audit.
No producer imports. Bits 0 and 1 say x/y beats each chain vertex; x->y.
"""
import itertools, json, pathlib, sys

def require(c, msg):
    if not c: raise RuntimeError(msg)

def graph(word):
    m=len(word); x=m; y=m+1
    a=[[False]*(m+2) for _ in range(m+2)]
    for i in range(m):
        for j in range(i+1,m): a[i][j]=True
        for v,k in ((x,0),(y,1)):
            a[v][i]=bool(word[i] & (1<<k)); a[i][v]=not a[v][i]
    a[x][y]=True
    return a

def literal(word):
    a=graph(word)
    return sum(all(a[p[i]][p[i+1]] for i in range(len(p)-1)) for p in itertools.permutations(range(len(a))))

def runs(word):
    total=0
    for A,B in ((0,1),(1,0)):
        for lab in itertools.product(range(3),repeat=len(word)):
            rows=[[i for i,c in enumerate(lab) if c==k] for k in range(3)]
            P,Q,R=rows
            if P and word[P[-1]]>>A&1: continue
            if Q:
                if not(word[Q[0]]>>A&1) or word[Q[-1]]>>B&1: continue
            elif A==1: continue
            if R and not(word[R[0]]>>B&1): continue
            total+=1
    return total

def subset(word):
    a=graph(word); n=len(a); d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if d[mask][i]:
                for j in range(n):
                    if not(mask>>j&1) and a[i][j]: d[mask|(1<<j)][j]+=d[mask][i]
    return sum(d[-1])

def main():
    checks=0; maxima=[]
    for m in range(5):
        maximum=0
        for w in itertools.product(range(4),repeat=m):
            h=literal(w); require(h==runs(w)==subset(w),f'count mismatch {w}')
            maximum=max(maximum,h); checks+=1
        maxima.append(maximum)
    # Deliberately erase the apex-orientation condition at the empty middle.
    require(runs(())==1,'empty source')
    require(literal((1,2))==5,'positive mixed control')
    out={'status':'PASS','all_word_lengths':list(range(5)),'word_cases':checks,'three_independent_counts':'literal permutations, generic adjacency subset DP, ordered three-run assignments','maxima':maxima,'empty_middle_hostile':'empty word has one x->y path, not two apex orders','scope':'All short words; does not claim global spectrum coverage.'}
    print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__': main()
