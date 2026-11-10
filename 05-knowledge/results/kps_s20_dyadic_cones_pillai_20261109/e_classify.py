# THM-5109 (E) cross-check over ABSTRACT states (no realisability assumed except h <= 2^j, 0 <= e <= h, e = R mod 2):
# enumerate every dyadic cone step through a transitive gadget (p >= 2, exit q < p) from states (h, e, 2^j) with
# j <= JMAX, using the closed forms of (B); then look for a second dyadic step from every output.
# Types: (i) q >= 2; (ii) q = 1.  A full sweep over all (h, e) is done for j <= JFULL; for JFULL < j <= JMAX only the
# states allowed by the classification (h in {1,5,13}, or the (ii-b) relation) are scanned, and the classification
# itself is re-verified on the full sweep.
import sys
JFULL=int(sys.argv[1]); JMAX=int(sys.argv[2]); MEXTRA=int(sys.argv[3])
ispow=lambda r: r>0 and r&(r-1)==0
def steps(h,e,j):
    R=1<<j; out=[]
    for p in range(2,j+MEXTRA+3):
        for q in range(1,p):
            m=p   # R' scales by 2^(m-p) or 2^(m-q); dyadic-ness of R' does not depend on m; e' depends on m only via p==m
            for mm in (p,p+1):
                hp=h+(1<<(p-1))*R
                ep=(1<<(q-1))*R if p==mm else 0
                if q>=2: Rp=(1<<(mm-q))*h+(1<<(mm-p+q))*3**(p-q-1)*R
                else:    Rp=(1<<(mm-1))*h-(1<<(mm-p))*e+(1<<(mm-p+1))*3**(p-2)*R
                if ispow(Rp): out.append(((mm,p,q),(hp,ep,Rp.bit_length()-1)))
    return out
def classified(h,e,j):
    if h in (1,5,13): return True
    for p,d in ((3,1),(5,5),(7,13)):
        if e==(1<<(p-1))*h-d*(1<<(j+1)): return True
    return False
one=[]; viol=0
for j in range(0,JFULL+1):
    for h in range(1,(1<<j)+1,2):
        for e in range(0,h+1):
            if (e-(1<<j))%2: continue
            s=steps(h,e,j)
            if s:
                one.append(((h,e,j),s))
                if not classified(h,e,j): viol+=1
print('full sweep j<=%d: states with a dyadic step: %d; classification violations: %d'%(JFULL,len(one),viol))
for st,s in one[:12]: print('  ',st,'->',sorted(set((g[1],g[2],o[0]) for g,o in s))[:6])
# extended scan j <= JMAX over classified states only
ext=[]
for j in range(JFULL+1,JMAX+1):
    cands=set()
    for h in (1,5,13):
        if h<=(1<<j):
            for e in range(0,h+1):
                if (e-(1<<j))%2==0: cands.add((h,e))
    for p,d in ((3,1),(5,5),(7,13)):
        for h in range(1,(1<<j)+1,2):
            e=(1<<(p-1))*h-d*(1<<(j+1))
            if 0<=e<=h and (e-(1<<j))%2==0: cands.add((h,e))
    for (h,e) in cands:
        s=steps(h,e,j)
        if s: ext.append(((h,e,j),s))
print('extended scan %d<j<=%d: states with a dyadic step: %d'%(JFULL,JMAX,len(ext)))
two=0; examples=[]
for st,s in one+ext:
    for g,o in s:
        if steps(*o):
            two+=1
            if len(examples)<5: examples.append((st,g,o,steps(*o)[:3]))
print('two-step chains (abstract):',two, examples)
