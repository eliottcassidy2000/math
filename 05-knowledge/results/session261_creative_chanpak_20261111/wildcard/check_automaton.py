#!/usr/bin/env python3
"""Exact positive automaton for tournaments with a specified transitive deletion."""
from itertools import permutations, product
from collections import deque
import random,json

def compile_automaton(k, arcs):
    states=[(pi,s) for pi in permutations(range(k)) for s in product(*([range(3)]*k+[range(2)]))]
    index={s:i for i,s in enumerate(states)}
    initial=[int(not any(s)) for pi,s in states]
    accept=[]
    transitions=[]
    for pi,s in states:
        ok=2 not in s
        for j in range(1,k):
            if not s[j] and not arcs[pi[j-1]][pi[j]]: ok=False
        accept.append(ok)
    for wordletter in range(1<<k):
        rows=[]
        for pi,s in states:
            destinations=[]
            for j in range(k+1):
                if not s[j] and j and not (wordletter>>pi[j-1]&1):continue
                ss=list(s)
                ss[j]=1 if j==k or not(wordletter>>pi[j]&1) else 2
                destinations.append(index[(pi,tuple(ss))])
            rows.append(destinations)
        transitions.append(rows)
    return initial,accept,transitions

def step(v,rows,cap=None):
    w=[0]*len(v)
    for a,dest in zip(v,rows):
        for j in dest:w[j]+=a
    if cap is not None:w=[min(x,cap) for x in w]
    return w

def automaton_count(k,arcs,word):
    v,a,tr=compile_automaton(k,arcs)
    for x in word:v=step(v,tr[x])
    return sum(x for x,b in zip(v,a) if b)

def tournament(k,arcs,word):
    n=k+len(word); adj=[[False]*n for _ in range(n)]
    for i in range(k):
        for j in range(k):adj[i][j]=arcs[i][j]
    for t,x in enumerate(word):
        for i in range(k):adj[i][k+t]=bool(x>>i&1);adj[k+t][i]=not adj[i][k+t]
        for s in range(t):adj[k+s][k+t]=True
    return adj

def hp(adj):
    n=len(adj)
    if not n:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1<<n):
        for i,x in enumerate(dp[mask]):
            if x:
                for j in range(n):
                    if not(mask>>j&1) and adj[i][j]:dp[mask|1<<j][j]+=x
    return sum(dp[-1])

def decide(k,arcs,target,budget=None):
    init,acc,tr=compile_automaton(k,arcs);cap=target+1
    init=tuple(min(x,cap) for x in init);q=deque([init]);prev={init:None}
    while q:
        v=q.popleft()
        if sum(x for x,b in zip(v,acc) if b)==target:
            w=[]
            while prev[v] is not None:v,c=prev[v];w.append(c)
            return {'status':'witness','word':w[::-1],'states':len(prev)}
        for c,rows in enumerate(tr):
            w=tuple(step(v,rows,cap))
            if w not in prev:prev[w]=(v,c);q.append(w)
        if budget and len(prev)>budget:return {'status':'budget','states':len(prev)}
    return {'status':'absent_all_lengths','states':len(prev)}

if __name__=='__main__':
    checks=0
    for k in range(4):
        pairs=list(__import__('itertools').combinations(range(k),2))
        for bits in range(1<<len(pairs)):
            arcs=[[False]*k for _ in range(k)]
            for b,(i,j) in enumerate(pairs):arcs[i][j]=bool(bits>>b&1);arcs[j][i]=not arcs[i][j]
            for n in range(7-k):
                words=list(product(range(1<<k),repeat=n))
                if len(words)>500:words=random.Random(261107+k+n).sample(words,500)
                for word in words:
                    a=automaton_count(k,arcs,word);b=hp(tournament(k,arcs,word))
                    assert a==b,(k,arcs,word,a,b)
                    checks+=1
    print(json.dumps({'independent_subset_dp_checks':checks}))
    for n in [1,3,5,7,9,13,21,29,31,43]:print(json.dumps({'k':1,'target':n,**decide(1,[[False]],n)}))
    for n in [7,21,29,43]:print(json.dumps({'k':2,'target':n,**decide(2,[[False,True],[False,False]],n,150000)}))
