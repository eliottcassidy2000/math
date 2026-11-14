# Dimension of the prefix-vector span in each phase of the interval-hub family (letters with y=0: 0,2; y=1: 1,3).
from ihub import T, fin, idx, NS
def fwd(v,let):
    w=[0]*NS
    for s in range(NS):
        if v[s]:
            for t in T[let][s]: w[t]+=v[s]
    return w
def rank(vecs,p=1000003):
    M=[[c%p for c in v] for v in vecs]; r=0
    for c in range(NS):
        pr=next((i for i in range(r,len(M)) if M[i][c]),None)
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]; inv=pow(M[r][c],p-2,p); M[r]=[(x*inv)%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(NS)]
        r+=1
    return r
v0=[0]*NS; v0[idx(0,0,0,0)]=1; v0[idx(1,0,0,0)]=1
pre=[v0]; fr=[v0]
for L in range(1,9):
    fr=[fwd(v,l) for v in fr for l in (0,2)]; pre+=fr
print("phase y=0 prefix span (from start):", rank(pre))
# then y=1 phase from those
mid=list(pre); fr=list(pre)
for L in range(1,6):
    fr=[fwd(v,l) for v in fr for l in (1,3)]; mid+=fr
print("after y=0 then y=1 phases:", rank(mid))
fin_all=list(mid); fr=list(mid[-2000:])
for L in range(1,5):
    fr=[fwd(v,l) for v in fr for l in (0,2)]; fin_all+=fr
print("after y=0, y=1, y=0 phases:", rank(fin_all))
# suffix (response) side: in y=0 phase from fin
def bwd(x,let): return [sum(x[t] for t in T[let][s]) for s in range(NS)]
suf=[fin[0],fin[1]]; fr=list(suf)
for L in range(1,9):
    fr=[bwd(x,l) for x in fr for l in (0,2)]; suf+=fr
print("suffix response span, y=0 phase:", rank(suf))
