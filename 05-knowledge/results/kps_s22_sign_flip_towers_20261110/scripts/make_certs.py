"""make_certs.py - one explicit certificate per odd c <= CMAX and per sign:
   a pair (B+, O+) with offset r = +1 (seed) or r = -1 (anti-seed) and slope M = c 2^v, of least order found.
Sources: (0) pairs (B, O) of order <= 9 from the census (portrait.c, files pt_m.txt);
         (1) one tower level on a census pair of order <= 9 of the opposite sign;
         (2) two tower levels on a census pair of order <= 6 of the same sign.
Line format:  sign c v2 order m code Omask levels      (levels = d,i,m;d,i,m or "-"), least order first.
usage: make_certs.py SIGN CMAX [c1,c2,...]   (optional list restricts the targets)"""
import sys
SIGN=int(sys.argv[1]); CMAX=int(sys.argv[2]); only=set(map(int,sys.argv[3].split(","))) if len(sys.argv)>3 else None
def odd(x):
    v=0
    while x%2==0: x//=2; v+=1
    return x,v
def step(s,l):
    r,g,M,E=s; d,i,m=l; x=i-1; al=2**d+1; ga=3**d+2**d
    return (-r, 2**(x+1+m+d)*M, al*2**m*g+2*ga*2**x*3**m*M-(E if x==0 else 0), al*2**x*M+((2*r-g) if m==0 else 0))
lv=[(d,i,m) for d in (1,2) for i in range(1,10) for m in range(0,5)]
best={}
def offer(c,v,order,prov,levels):
    if c>CMAX or (only is not None and c not in only): return
    k=(order,v)
    if c not in best or k<best[c][0]: best[c]=(k,prov,levels)
for mm in range(1,10):
    for l in open("pt_%d.txt"%mm):
        m_,code,O,r,h,M,E=map(int,l.split())
        if M<=0: continue
        s=(r,h+r,M,E); prov=(m_,code,O)
        if r==SIGN:
            c,v=odd(M); offer(c,v,m_,prov,())
            if m_<=6:
                for l1 in lv:
                    s1=step(s,l1)
                    if s1[2]<=0: continue
                    for l2 in lv:
                        s2=step(s1,l2)
                        if s2[2]>0:
                            c,v=odd(s2[2]); offer(c,v,m_+sum(z[1]+z[0]+z[2]+1 for z in (l1,l2)),prov,(l1,l2))
        elif m_<=9:
            for l1 in lv:
                s1=step(s,l1)
                if s1[2]>0:
                    c,v=odd(s1[2]); offer(c,v,m_+l1[1]+l1[0]+l1[2]+1,prov,(l1,))
targets=[c for c in range(1,CMAX+1,2) if only is None or c in only]
miss=[c for c in targets if c not in best]
for c in targets:
    if c in best:
        (order,v),prov,levels=best[c]
        print(SIGN,c,v,order,prov[0],prov[1],prov[2],";".join("%d,%d,%d"%z for z in levels) if levels else "-")
print("# sign %d: certificates for %d of %d odd c <= %d ; missing: %s"%(SIGN,len(targets)-len(miss),len(targets),CMAX,miss),file=sys.stderr)
