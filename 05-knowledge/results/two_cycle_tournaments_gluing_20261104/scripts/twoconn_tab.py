import sys
from collections import Counter
def rows_from_d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    rows = [0] * n
    for i in range(n):
        for j in range(n):
            if bits[i * n + j]: rows[i] |= 1 << j
    return rows
def strong(rows, verts):
    vs=list(verts); s0=vs[0]
    for fwd in (True,False):
        seen={s0}; fr=[s0]
        while fr:
            nf=[]
            for v in fr:
                for u in vs:
                    if u not in seen and ((rows[v]>>u)&1 if fwd else (rows[u]>>v)&1): seen.add(u); nf.append(u)
            fr=nf
        if len(seen)<len(vs): return False
    return True
K=int(sys.argv[1])
for n in range(int(sys.argv[2]), int(sys.argv[3])+1):
    c=Counter(); t=Counter()
    for l in open('clsK%d_%d.txt'%(K,n)):
        h,d6=l.split(); h=int(h); rows=rows_from_d6(d6); t[h]+=1
        if all(strong(rows,[v for v in range(n) if v!=x]) for x in range(n)): c[h]+=1
    print("order %d (K=%d): classes by hc %s; 2-connected by hc %s"%(n,K,dict(sorted(t.items())),dict(sorted(c.items()))), flush=True)
