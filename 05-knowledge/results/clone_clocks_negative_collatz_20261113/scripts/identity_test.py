"""Test: for a vertex v with matching number nu_v = 1, is X_1 - 2 X_2 = X_0, where X_i = number of Hamiltonian
words visiting v exactly i times (v a module position), X_0 = H(S - v)?   (single module v; no u)"""
import random, sys
sys.path.insert(0,'scripts')
from funnel import H_only
def words_by_visits(adj, n, v, cap=4):
    others=[x for x in range(n) if x!=v]; pos={x:i for i,x in enumerate(others)}; full=(1<<len(others))-1
    from collections import defaultdict
    fr=defaultdict(int)
    fr[(0,1,v)]+=1
    for x in others: fr[(1<<pos[x],0,x)]+=1
    tot=defaultdict(int)
    while fr:
        nf=defaultdict(int)
        for (m,c,last),val in fr.items():
            if m==full: tot[c]+=val
            for w in range(n):
                if not (adj[last]>>w&1): continue
                if w==v:
                    if last==v or c==cap: continue
                    nf[(m,c+1,v)]+=val
                else:
                    if m>>pos[w]&1: continue
                    nf[(m|(1<<pos[w]),c,w)]+=val
        fr=nf
    return dict(tot)
def matching(adj,n,v):
    O=[x for x in range(n) if x!=v and adj[v]>>x&1]; I=[x for x in range(n) if x!=v and not adj[v]>>x&1]
    # max matching size via brute force augmenting
    match={}
    def try_o(o,seen):
        for i in I:
            if adj[o]>>i&1 and i not in seen:
                seen.add(i)
                if i not in match or try_o(match[i],seen):
                    match[i]=o; return True
        return False
    return sum(1 for o in O if try_o(o,set()))
if __name__ == '__main__':
    rng=random.Random(5)
    stats={}
    for trial in range(3000):
        n=rng.randint(5,8)
        adj=[0]*n
        for i in range(n):
            for j in range(i+1,n):
                if rng.random()<0.5: adj[i]|=1<<j
                else: adj[j]|=1<<i
        for v in range(n):
            nu=matching(adj,n,v)
            if nu!=1: continue
            t=words_by_visits(adj,n,v)
            X0=t.get(0,0); X1=t.get(1,0); X2=t.get(2,0)
            # X0 as computed = words with 0 v-symbols = H(S-v)
            key=(X1-2*X2==X0)
            stats[key]=stats.get(key,0)+1
            if not key and stats[key]<=5: print('counterexample n',n,'v',v,'X0,X1,X2',X0,X1,X2, 'X1-2X2-X0=',X1-2*X2-X0)
    print(stats)
