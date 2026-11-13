import sys; sys.path.insert(0,'scripts')
from funnel import H_only
def transitive_with_back(order_len, back):
    """vertices 0..n-1 in transitive order (i->j for i<j) except arcs in back (pairs (hi,lo)) reversed to hi->lo"""
    n = order_len; adj=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if (j,i) in back: adj[j]|=1<<i
            else: adj[i]|=1<<j
    return adj
def blow(adj,n,u,m):
    from tern_check import blowup
    A,n2,block,newid = blowup(adj,n,u,m); return A,n2
def fit3(vals):
    # vals[j-1] for j=1..5 ; fit a*3^j + b*2^j + c
    import sympy as sp
    a,b,c=sp.symbols('a b c')
    sol=sp.solve([a*3**j+b*2**j+c-vals[j-1] for j in (1,2,3)],[a,b,c])
    ok=all(sol[a]*3**j+sol[b]*2**j+sol[c]==vals[j-1] for j in range(1,len(vals)+1))
    return (sol[a],sol[b],sol[c],ok)
for pre in (0,1,2):          # extra transitive vertices between the two anchors (like b)
  for p in range(1,4):
    for q in range(1,4):
        # order: a, [pre vertices], c, u, D(p), F(q)
        names=['a']+['b%d'%i for i in range(pre)]+['c','u']+['d%d'%i for i in range(p)]+['f%d'%i for i in range(q)]
        idx={nm:i for i,nm in enumerate(names)}
        back=set()
        for i in range(p): back.add((idx['d%d'%i], idx['a']))
        for i in range(q): back.add((idx['f%d'%i], idx['c']))
        adj=transitive_with_back(len(names), back)
        vals=[H_only(*blow(adj,len(names),idx['u'],m)) for m in range(1,5)]
        print('pre',pre,'p',p,'q',q, 'H(j)=',vals, ' fit 3^j,2^j,1:', fit3(vals))
