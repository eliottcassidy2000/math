"""Two-module visit numbers a[i][i'] (u visited i times, v visited i' times) by literal word DP, and the
two-clock test  H(j,k) = 3^j (alpha 2^k + beta) - 1."""
import sys; sys.path.insert(0,'scripts')
from functools import lru_cache
def visit_table(adj, n, u, v, CU=6, CV=6):
    others = [x for x in range(n) if x not in (u, v)]
    pos = {x:i for i,x in enumerate(others)}
    full = (1 << len(others)) - 1
    # dp over (mask, cu, cv, last) ; last in others or 'u' or 'v'
    from collections import defaultdict
    dp = defaultdict(int)
    for x in range(n):
        if x == u: dp[(0,1,0,u)] += 1
        elif x == v: dp[(0,0,1,v)] += 1
        else: dp[(1<<pos[x],0,0,x)] += 1
    # process in order of total length
    by_len = defaultdict(list)
    def L(state): m,cu,cv,last = state; return bin(m).count('1')+cu+cv
    frontier = dict(dp)
    total = defaultdict(int)
    while frontier:
        nf = defaultdict(int)
        for (m,cu,cv,last),val in frontier.items():
            if m == full: total[(cu,cv)] += val
            for w in range(n):
                if not (adj[last] >> w & 1): continue
                if w == u:
                    if last == u or cu == CU: continue
                    nf[(m,cu+1,cv,u)] += val
                elif w == v:
                    if last == v or cv == CV: continue
                    nf[(m,cu,cv+1,v)] += val
                else:
                    if m >> pos[w] & 1: continue
                    nf[(m | (1<<pos[w]),cu,cv,w)] += val
        frontier = nf
    return dict(total)
def twoclock_params(tab):
    """return (alpha, beta) if the table is a two-clock certificate, else None"""
    g = lambda i,ip: tab.get((i,ip),0)
    for (i,ip),val in tab.items():
        if val and (i >= 4 or ip >= 3 or i == 0 or ip == 0): 
            if i==0 or ip==0: continue
            return None
    A = g(3,1); al = g(3,2)
    if A <= 0 or al <= 0: return None
    if (g(2,1),g(1,1)) != (3*A, 3*A-1): return None
    if (g(2,2),g(1,2)) != (3*al, 3*al): return None
    return (al, A - 2*al)
if __name__ == '__main__':
    from tern_check import parse
    import itertools
    sk = [('1110011111111111001111111111',8,3),('1110001111111111111100111111',8,6)]
    for s,n,u in sk:
        adj = parse(s,n)
        for v in range(n):
            if v == u: continue
            tab = visit_table(adj,n,u,v)
            print(s,'u',u,'v',v, {k:tab[k] for k in sorted(tab) if tab[k]}, '->', twoclock_params(tab))
