import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import itertools, subprocess, sys
from collections import Counter
n=int(sys.argv[1])
lines=subprocess.run(['gentourng',str(n)],capture_output=True,text=True).stdout.split()
def parse(s):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if s[k]=='1': out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def paths(out):
    res=[]
    def dfs(p,used):
        if len(p)==n: res.append(tuple(p)); return
        nb=out[p[-1]]&~used
        while nb:
            u=(nb&-nb).bit_length()-1; nb&=nb-1
            dfs(p+[u],used|1<<u)
    for s in range(n): dfs([s],1<<s)
    return res
def allpaths(out):
    # all nonempty paths as tuples
    res=[]
    def dfs(p,used):
        res.append(tuple(p))
        nb=out[p[-1]]&~used
        while nb:
            u=(nb&-nb).bit_length()-1; nb&=nb-1
            dfs(p+[u],used|1<<u)
    for s in range(n): dfs([s],1<<s)
    return res
stats=Counter(); ex=[]; mE={}
for s in lines:
    if s.startswith('>'): continue
    out=parse(s)
    HP=paths(out); H=len(HP)
    AP=allpaths(out)
    bymask={}
    for p in AP:
        m=0
        for v in p: m|=1<<v
        bymask.setdefault(m,[]).append(p)
    full=(1<<n)-1
    for X in range(1,full):
        inX=lambda v:(X>>v)&1
        Z=E=0
        for p in HP:
            t=sum(1 for i in range(n-1) if inX(p[i]) and not inX(p[i+1]))
            if t==0: Z+=1
            else: E+=t-1
        R=0
        for m,ps in bymask.items():
            if m==full: continue
            cm=full^m
            for q1 in ps:
                u=q1[-1]
                if inX(u): continue
                for q2 in bymask.get(cm,[]):
                    w=q2[0]
                    if inX(w) and (out[w]>>u)&1: R+=1
        d=E+R-Z
        if E>=1:
            stats[d]+=1
            if d==3: ex.append((H,E,R,Z))
            mE[E]=min(mE.get(E,99),d)
print('n',n,'delta distribution over double-crossing cuts:',sorted(stats.items())[:10])
print('extremal (H,E,R,Z) counts:',Counter(ex)); print('min delta by E:',sorted(mE.items())[:8])
# structure of extremal cases: (E,R,Z) histogram over delta==3 with E>=1, and min delta per E
