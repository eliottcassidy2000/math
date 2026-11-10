# P4: exhaustive two-step transitive dyadic cone chains from every dyadic base of order <= 8.
# Transitive gadget TT_m, entry position p >= 2, exit position q < p (exits at or after the entry give 2^(m-q) h',
# never a power of two because h' > 1).  Closed forms (THM-5109 (A)):
#   h' = h + 2^(p-1) R
#   e' = 2^(q-1) R if p == m else 0            (q < p)
#   R' = 2^(m-q) h + 2^(m-p+q) 3^(p-q-1) R     (2 <= q < p)
#   R' = 2^(m-1) h - 2^(m-p) e + 2^(m-p+1) 3^(p-2) R   (q = 1 < p)
import sys
PORTS=sys.argv[1]; MMAX=int(sys.argv[2])
ispow=lambda r: r>0 and r&(r-1)==0
bases=set()
for line in open(PORTS):
    p=list(map(int,line.split()))
    n,code,u,h,E,O,hu=p[:7]
    R=E+O; e=O-E
    if ispow(R): bases.add((h,e,R))
def steps(h,e,R):
    out=[]
    for m in range(2,MMAX+1):
        for p in range(2,m+1):
            hp=h+(1<<(p-1))*R
            for q in range(1,p):
                ep=(1<<(q-1))*R if p==m else 0
                if q>=2: Rp=(1<<(m-q))*h+(1<<(m-p+q))*3**(p-q-1)*R
                else:    Rp=(1<<(m-1))*h-(1<<(m-p))*e+(1<<(m-p+1))*3**(p-2)*R
                if ispow(Rp): out.append(((m,p,q),(hp,ep,Rp)))
    return out
one=0; two=[]
firststeps={}
for b in sorted(bases):
    s1=steps(*b)
    if s1: one+=1; firststeps[b]=s1
    for g1,st in s1:
        for g2,st2 in steps(*st):
            two.append((b,g1,st,g2,st2))
print('dyadic bases (distinct triples, order<=8):',len(bases))
print('bases with a dyadic one-step transitive cone (m<=%d):'%MMAX,one)
for b,s in sorted(firststeps.items())[:40]:
    print('  base',b,'->',[(g,st[0]) for g,st in s][:6])
print('two-step chains:',len(two))
for c in two[:20]: print('  ',c)
