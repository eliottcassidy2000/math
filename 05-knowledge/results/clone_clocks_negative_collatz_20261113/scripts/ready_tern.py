"""For digit-1 ternary clocks (S,u) from the back-arc scans: does any vertex a have C2 digit 1 (either birth
convention) in S[u<-TT_j] for j = 1,2,3?  (then the clock is READY at every j, by the exponential-polynomial
structure once the support bound is in place; here we only screen)."""
import sys, glob, re
sys.path.insert(0,'scripts')
from digits import digits
from tern_check import blowup
def build(n, back):
    adj=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if (j,i) in back: adj[j]|=1<<i
            else: adj[i]|=1<<j
    return adj
seen=set(); res=[]
files = sys.argv[1:]
for f in files:
    for line in open(f):
        if not line.startswith('TERN') or ' sigma=1 ' not in line: continue
        n=int(re.search(r'n=(\d+)',line).group(1)); u=int(re.search(r' u=(\d+)',line).group(1))
        A=int(re.search(r' A=(\d+)',line).group(1))
        back=tuple(sorted(tuple(map(int,p.split('>'))) for p in re.search(r'back=([\d>,]+)',line).group(1).strip(',').split(',')))
        key=(n,u,back)
        if key in seen: continue
        seen.add(key)
        adj=build(n,set(back))
        H,sig,tau=digits(adj,n)
        assert H==3*A-1, (H,A)
        cand=[('sig',a) for a in range(n) if sig[a]==1]+[('tau',a) for a in range(n) if tau[a]==1]
        if not cand: continue
        ok=[]
        for (conv,a) in cand:
            good=True
            for m in (2,3):
                A2,n2,block,newid=blowup(adj,n,u,m)
                H2,s2,t2=digits(A2,n2)
                aa=newid[a] if a!=u else None
                if aa is None: good=False; break
                d=(s2 if conv=='sig' else t2)[aa]
                if d!=1: good=False; break
            if good: ok.append((conv,a))
        res.append((A,n,u,back,cand,ok))
        print('A',A,'n',n,'u',u,'back',back,'digit-1 at j=1:',cand,' persistent:',ok, flush=True)
print('checked',len(seen),'clocks; with a digit-1 vertex at j=1:',len(res),'; persistent:',sum(1 for r in res if r[5]))
