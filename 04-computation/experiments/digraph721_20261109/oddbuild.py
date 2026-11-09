"""Explicit odd realisers on exactly L(k) vertices for every odd k other than 1, 5, 37, following the THM-5056
induction; brute-force check of d and of an adjacent-twin pair (a certificate of oddness, THM-4996 B2.1)."""
import sys, random
CENSUS='dom_K511.out'
def L(k): return k.bit_length()
def dcount(adj):
    n=len(adj); cl=[adj[v]|(1<<v) for v in range(n)]; full=(1<<n)-1; c=0
    dom=[0]*(1<<n)
    for S in range(1,1<<n):
        b=(S&-S).bit_length()-1; dom[S]=dom[S&(S-1)]|cl[b]
        if dom[S]==full: c+=1
    return c
def pad(adj,k=1): return adj+[0]*k
def padto(adj,n): return adj+[0]*(n-len(adj))
def cone(adj): n=len(adj); return [a|(1<<n) for a in adj]+[(1<<n)-1]
def union(a,b): na=len(a); return a+[x<<na for x in b]
def join(a,b):
    na,nb=len(a),len(b); fa=(1<<na)-1; fb=((1<<nb)-1)<<na
    return [x|fb for x in a]+[(y<<na)|fa for y in b]
def K(n): return [((1<<n)-1)^(1<<v) for v in range(n)]
def ph(adj,u,l):
    n=len(adj); h=n; new=list(adj); hm=0
    for v in range(n):
        if v!=u: new[v]|=1<<h; hm|=1<<v
    for i in range(l): hm|=1<<(n+1+i)
    return new+[hm]+[1<<h]*l
def threshold(k):
    adj=[0]
    for i in range(1,L(k)):
        adj = cone(adj) if (k>>i)&1 else pad(adj)
    return adj
def twins(adj):
    n=len(adj)
    for x in range(n):
        for y in range(x+1,n):
            if (adj[x]>>y)&1 and (adj[x]|(1<<x))==(adj[y]|(1<<y)): return (x,y)
    return None
# base: least-order odd realisers below 512 from the census (connected odd class, possibly times a threshold graph)
conn={}
for line in open(CENSUS):
    p=line.split(); n=int(p[0]); d=int(p[2])
    if p[3]!='O': continue
    if d in conn and len(conn[d])<=n: continue
    e=p[5] if len(p)>5 else ''; adj=[0]*n
    for i in range(0,len(e),2):
        u=int(e[i],16); v=int(e[i+1],16); adj[u]|=1<<v; adj[v]|=1<<u
    conn[d]=adj
base={}
for k in range(3,512,2):
    best=None
    for a,g in conn.items():
        if k%a: continue
        m=k//a; G=g if m==1 else union(g,threshold(m))
        if best is None or len(G)<len(best): best=G
    if best is not None: base[k]=best
P3=[0b110,0b001,0b001]; E2=[0,0]
K4e=[0b1110,0b1101,0b0011,0b0011]; C4=[0b1010,0b0101,0b1010,0b0101]
SEEDS={5:[(K4e,2)], 1:[(union(K(2),C4),2)], 37:[(join(K(2),union(P3,E2)),5),(union(K(2),threshold(55)),None)]}
def deg(adj): return [bin(a).count('1') for a in adj]
def first_ph(G,w):
    if w is None: w=[v for v in range(len(G)) if deg(G)[v]==2 and v>1][0]
    return ph(G,w,2)
def offset_graph(c,n):
    """good odd graph on n vertices with d = 2^(n-1) + c, built from a seed by PH(3,2) then moves +2 (R) / +3 (PH(3,2))."""
    for G0,w in SEEDS[c]:
        m=len(G0)+3
        if n<m: continue
        r=n-m
        # write r = 2a + 3b with a,b >= 0 (r != 1)
        if r==1: continue
        b=r%2*1; a=(r-3*b)//2
        G=first_ph(G0,w)
        tw=twins(G)
        for _ in range(a):
            dg=deg(G); p=[v for v in range(len(G)) if dg[v]==1 and v not in tw][0]; G=ph(G,p,1)
        for _ in range(b):
            dg=deg(G); w2=[v for v in range(len(G)) if dg[v]==2 and v not in tw][0]; G=ph(G,w2,2)
        return G
    return None
memo={}
def build(k):
    if k in memo: return memo[k]
    Lk=L(k)
    if k<512: G=base[k]
    elif (k>>(Lk-2))&1:                        # leading 11: K_2 v G
        G=join(K(2),padto(threshold(k-3*2**(Lk-2)),Lk-2)) if k-3*2**(Lk-2)>0 else None
    else:
        x=k-2**(Lk-1)
        if x not in (1,5,37): G=cone(padto(build(x),Lk-1))
        elif k==513: G=union(K(2),threshold(171))
        else: G=offset_graph(x,Lk)
    assert G is not None and len(G)==Lk, (k,None if G is None else len(G))
    memo[k]=G; return G
lo,hi=int(sys.argv[1]),int(sys.argv[2]); samples=int(sys.argv[3]) if len(sys.argv)>3 else 0
ks=list(range(lo|1,hi,2))
if samples: random.seed(5); ks=sorted(random.sample(ks,samples))
out=open(sys.argv[4],'w') if len(sys.argv)>4 else None
bad=0; cnt=0
for k in ks:
    if k in (1,5,37): continue
    G=build(k); cnt+=1
    if len(G)!=L(k): bad+=1; print('BAD order',k)
    if out:
        n=len(G); edges=' '.join(f'{u}-{v}' for u in range(n) for v in range(u+1,n) if (G[u]>>v)&1)
        out.write(f'odd{k} O {n} {k} | {edges}\n')
print('built',cnt,'values in [%d,%d) on exactly L(k) vertices; order failures=%d'%(lo,hi,bad))
