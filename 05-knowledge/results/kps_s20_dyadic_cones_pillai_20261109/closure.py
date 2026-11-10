# dyadic closure: triples (h,e,R) with R a power of two, closed under one-arc gadget matrices (exit in gadget)
# and S-padding (h,e,R)->(h,0,2R). Only outputs whose R is a power of two are kept.
import sys
GMAX=int(sys.argv[1]); HMAX=int(sys.argv[2]); RMAXF=int(sys.argv[3]) if len(sys.argv)>3 else 64
mats=set()
for line in open('gadgets5.txt'):
    p=list(map(int,line.split()))
    if p[0]<=GMAX: mats.add(tuple(p[4:]))
mats=sorted(mats)
ispow=lambda r: r>0 and r&(r-1)==0
start=[(1,1,1)]
seen=set(); frontier=[]
def add(s):
    h,e,R=s
    if h>HMAX or R>RMAXF*HMAX: return
    if s not in seen: seen.add(s); frontier.append(s)
add((1,1,1))
while frontier:
    s=frontier.pop()
    h,e,R=s
    add((h,0,2*R))
    for m in mats:
        h2=m[0]*h+m[1]*e+m[2]*R; e2=m[3]*h+m[4]*e+m[5]*R; R2=m[6]*h+m[7]*e+m[8]*R
        if ispow(R2): add((h2,e2,R2))
vals=sorted(set(h for h,e,R in seen))
print('gadget order <=',GMAX,' distinct matrices',len(mats),' dyadic triples',len(seen),' values',len(vals))
print('values:',vals[:80])
