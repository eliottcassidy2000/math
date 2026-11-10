import subprocess, itertools
EXE='./resp.exe'
def run(tours):
    """tours: list of out-masks lists. returns list of (H, Gamma, [R...])"""
    inp=''.join('%d %s\n'%(len(t),' '.join(map(str,t))) for t in tours)
    out=subprocess.run([EXE],input=inp,capture_output=True,text=True).stdout.split('\n')
    res=[]
    for line in out:
        if not line.strip(): continue
        p=list(map(int,line.split('|')[0].split())); res.append((p[0],p[1],p[2:]))
    return res
def from_arcs(n,arcs):
    out=[0]*n
    for a,b in arcs: out[a]|=1<<b
    return out
def transitive_with(n,back):
    """vertices 0..n-1, i->j for i<j except pairs in back (j->i)"""
    out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if (i,j) in back or (j,i) in back: out[j]|=1<<i
            else: out[i]|=1<<j
    return out
def add_vertex(out,beats):
    """new vertex beats set 'beats', beaten by all others"""
    n=len(out); new=list(out)
    m=0
    for v in range(n):
        if v in beats: m|=1<<v
        else: new[v]|=1<<n
    new.append(m); return new
def add_source(out):
    n=len(out); return list(out)+[(1<<n)-1]
def ispow(r): return r>0 and r&(r-1)==0
def nearjoin(out,u,m):
    """B-hat at port u (w beats only u), then m+1 vertices each beating only w. value h + 2^(m+1) R_u"""
    t=add_vertex(out,{u}); w=len(out)
    for i in range(m+1): t=add_vertex(t,{w})
    return t

def run2(tours):
    """returns list of (H, Gamma, [R...], [e...])"""
    inp=''.join('%d %s\n'%(len(t),' '.join(map(str,t))) for t in tours)
    out=subprocess.run([EXE],input=inp,capture_output=True,text=True).stdout.split('\n')
    res=[]
    for line in out:
        if not line.strip(): continue
        a,b=line.split('|'); p=list(map(int,a.split())); ee=list(map(int,b.split()))
        res.append((p[0],p[1],p[2:],ee))
    return res
def join(T,v,G,u):
    """N = T => G with u -> v reversed. T,G out-mask lists. returns out-masks; G vertices offset by len(T)"""
    n=len(T); g=len(G); out=[0]*(n+g)
    for i in range(n): out[i]=T[i] | (((1<<g)-1)<<n)
    for j in range(g): out[n+j]=G[j]<<n
    out[v]&=~(1<<(n+u)); out[n+u]|=1<<v
    return out
def TT(m): return [sum(1<<k for k in range(i+1,m)) for i in range(m)]
