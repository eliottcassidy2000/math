import sys; sys.path.insert(0,'scripts')
from funnel import hp_stats, is_tournament
def parse(s, n):
    adj=[0]*n; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if s[idx]=='1': adj[i]|=1<<j
            else: adj[j]|=1<<i
            idx+=1
    return adj
def blowup(adj, n, u, m):
    """replace vertex u by TT_m (u_1 -> u_2 -> ... -> u_m), uniform external arcs. returns adj2, n2, list of block vertices"""
    others=[v for v in range(n) if v!=u]
    newid={v:i for i,v in enumerate(others)}
    n2=len(others)+m
    block=list(range(len(others), n2))
    A=[0]*n2
    for v in others:
        for w in others:
            if adj[v]>>w&1: A[newid[v]]|=1<<newid[w]
        for b in block:
            if adj[v]>>u&1: A[newid[v]]|=1<<b
            else: A[b]|=1<<newid[v]
    for i,b in enumerate(block):
        for b2 in block[i+1:]: A[b]|=1<<b2
    return A, n2, block, newid
def describe(s, n, u):
    adj=parse(s,n)
    print('S =', s, 'u =', u)
    for v in range(n):
        print('  ', v, 'out:', [w for w in range(n) if adj[v]>>w&1], 'outdeg', bin(adj[v]).count('1'))
    for m in range(1,6):
        A,n2,block,newid=blowup(adj,n,u,m)
        H,s_,e_=hp_stats(A,n2)
        # readiness candidates: vertices y with indegree 1 (sole predecessor x) and delta_y = H(T-y) - s_y = s_x(T-y)
        ready=[]
        for y in range(n2):
            preds=[x for x in range(n2) if A[x]>>y&1]
            if len(preds)==1:
                x=preds[0]
                # build T-y
                keep=[v for v in range(n2) if v!=y]; idm={v:i for i,v in enumerate(keep)}
                B=[0]*(n2-1)
                for v in keep:
                    for w in keep:
                        if A[v]>>w&1: B[idm[v]]|=1<<idm[w]
                Hy,sy,ey=hp_stats(B,n2-1)
                ready.append((y, x, sy[idm[x]]))
        print('  j=%d  H=%d  (H+1)/3^j=%s  indegree-1 roots (y,pred,delta):'%(m,H,(H+1)/3**m), ready)
if __name__=='__main__':
    for s,n,u in [('1110011111111111001111111111',8,3),('111001111111111111001111111111111110',9,3),('111010111111111111100111111111111111',9,3),('111100111111111111100111111111111111',9,3)]:
        describe(s,n,u)
