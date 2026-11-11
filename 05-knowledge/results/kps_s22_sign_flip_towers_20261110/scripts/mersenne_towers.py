"""mersenne_towers.py - depth-3 pure towers with offset -1 and slope a power of two (rays 2^s - 1), found by the state
recursion in the E1 box and verified by the transfer-matrix count of the explicit tournaments (kind-pasteur-2026-11-10-S22)."""
from tower_build import build, adjacency, count_apex, start_state, step
XMAX,MMAX,KMAX=10,4,10
starts=[(k+1,[1]) for k in range(KMAX+1)]+[(1,[])]
lv=[(d,x+1,m) for d in (1,2) for x in range(XMAX+1) for m in range(MMAX+1)]
hits=[]
for st in starts:
    s0=start_state(st)
    for l1 in lv:
        s1=step(s0,l1)
        for l2 in lv:
            s2=step(s1,l2)
            for l3 in lv:
                r,M,h,E=step(s2,l3)
                if M>0 and M&(M-1)==0: hits.append((st[0]+sum(i+d+m+1 for d,i,m in (l1,l2,l3))+2,M.bit_length()-1,st,(l1,l2,l3)))
hits.sort(key=lambda z:(z[0],z[1],z[2][0]))
print("depth-3 states with r = -1 and M a power of two in the E1 box:",len(hits))
ok=True
for order,v,st,levels in hits:
    line="  order(t=1) %d  M = 2^%d  start TT_%d, out-set %s  levels (delta,i,m) = %s"%(order,v,st[0],st[1],levels)
    for t in (1,2):
        N,q,beats=build(st,list(levels),t); A=adjacency(N,q,beats); Hc=count_apex(N,q,A); ok&=(Hc==2**(v+t)-1 and len(A)==order+t-1)
        line+="  | t=%d: order %d, H = %d = 2^%d - 1: %s"%(t,len(A),Hc,v+t,Hc==2**(v+t)-1)
    print(line)
print("ALL VERIFIED" if ok and hits else "MISMATCH")
