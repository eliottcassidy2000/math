# Independent exact check (rational arithmetic on class representatives, explicit paths):
# for every unit class c mod 3^7 with c mod 27 not in {1,14}, find a legal reverse E-path of length<=6
# (legality decided by c mod 3^7) whose multiplier 2^K/3^s <= 8/9 at some prefix.  Print the certificate table.
from fractions import Fraction as F
R=6; M=3**(R+1)
def legal_paths(c, depth):
    # yields (path list of k, K, s) for all legal paths up to depth, using representative c (valid mod 3^(R+1-s))
    stack=[(c, 0, 0, [])]
    while stack:
        x,K,s,path=stack.pop()
        if s>0: yield path,K,s
        if s==depth: continue
        prec=R+1-s
        if prec<2: continue
        mod=3**prec
        xr=x%mod
        k0=0 if xr%3==1 else 1
        for k in range(k0, 2*depth+4, 2):
            y=(2**k*x-1)
            if y%3: raise SystemExit("bad")
            y//=3
            if y%3==0: continue
            stack.append((y,K+k,s+1,path+[k]))
worst=F(0); worstc=None; cert={}
for c in range(M):
    if c%3==0 or c%27 in (1,14): continue
    best=None
    for path,K,s in legal_paths(c+M, R):   # representative c+M (>1) to avoid degenerate small values
        r=F(2**K,3**s)
        if best is None or r<best[0]: best=(r,path)
    assert best is not None and best[0]<=F(8,9), (c,best)
    cert[c]=best
    if best[0]>worst: worst=best[0]; worstc=c
print("classes checked:",len(cert)," worst best-multiplier:",worst,"at class",worstc,"path",cert[worstc][1])
from collections import Counter
print("multiplier histogram (best per class):",sorted(Counter(str(v[0]) for v in cert.values()).items(),key=lambda t:F(t[0]))[:12])
