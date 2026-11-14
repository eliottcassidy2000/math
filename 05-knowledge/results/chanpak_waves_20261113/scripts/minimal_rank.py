import itertools, random
NS=36
def idx(o,a,b,d): return ((o*3+a)*3+b)*2+d
T=[[None]*NS for _ in range(4)]
for let in range(4):
    xi,yi=let>>1,let&1
    for o in range(2):
        for a in range(3):
            for b in range(3):
                for d in range(2):
                    s=idx(o,a,b,d); h1=yi if o else xi; h2=xi if o else yi
                    succ=[idx(o,2 if h1 else 1,b,d)]
                    if b or h1: succ.append(idx(o,a,2 if h2 else 1,d))
                    if d or h2: succ.append(idx(o,a,b,1))
                    T[let][s]=succ
fin=[[0]*NS for _ in range(2)]
for e in range(2):
    for o in range(2):
        for a in range(3):
            for b in range(3):
                for d in range(2):
                    h12=(1-e) if o else e
                    fin[e][idx(o,a,b,d)]=0 if (a==2 or b==2 or (b==0 and not h12)) else 1
def fwd(v,let):
    w=[0]*NS
    for s in range(NS):
        if v[s]:
            for t in T[let][s]: w[t]+=v[s]
    return w
def bwd(x,let):  # response: x'[s]=sum_{t in T[let][s]} x[t]
    return [sum(x[t] for t in T[let][s]) for s in range(NS)]
def rank_mod(vecs,p):
    rows=[[c%p for c in v] for v in vecs]; r=0; piv_cols=[]
    M=[row[:] for row in rows]; ncol=NS
    for c in range(ncol):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]%p: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]; inv=pow(M[r][c],p-2,p) if p>2 else 1
        M[r]=[(x*inv)%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(ncol)]
        r+=1
    return r
v0=[0]*NS; v0[idx(0,0,0,0)]=1; v0[idx(1,0,0,0)]=1
pref=[v0]; frontier=[v0]
for L in range(1,7):
    frontier=[fwd(v,l) for v in frontier for l in range(4)]; pref+=frontier
suf=[fin[0],fin[1]]; fr=[fin[0],fin[1]]
for L in range(1,7):
    fr=[bwd(x,l) for x in fr for l in range(4)]; suf+=fr
for p in (3,5,7,1000003):
    print("p=%d: prefix span rank %d, suffix span rank %d"%(p,rank_mod(pref,p),rank_mod(suf,p)))
# Hankel rank: matrix H[i][j]=pref_i . suf_j over a sample
import random
random.seed(1)
P=random.sample(pref,min(400,len(pref))); S=random.sample(suf,min(400,len(suf)))
Hk=[[sum(a*b for a,b in zip(pv,sv)) for sv in S] for pv in P]
def rank_rows(M,p):
    M=[[c%p for c in row] for row in M]; r=0; ncol=len(M[0])
    for c in range(ncol):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]; inv=pow(M[r][c],p-2,p)
        M[r]=[(x*inv)%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(ncol)]
        r+=1
    return r
for p in (3,5,1000003):
    print("Hankel rank (sample 400x400) mod %d: %d"%(p,rank_rows(Hk,p)))
