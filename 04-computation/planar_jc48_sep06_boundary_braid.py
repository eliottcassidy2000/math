#!/usr/bin/env python3
"""Rational Rouche certificate for actual quartet braids; candidate under audit.

Default: verify a frozen rational path certificate without floating arithmetic.
--produce: construct candidate paths using NumPy, accepting only exact tubes.
"""
from fractions import Fraction as Q
from math import comb
from itertools import permutations, combinations, product
from pathlib import Path
import json
import sys
import hashlib
import gzip


class C:
    __slots__=('r','i')
    def __init__(self,r=0,i=0): self.r,self.i=Q(r),Q(i)
    def __add__(self,b):
        if not isinstance(b,C): b=C(b)
        return C(self.r+b.r,self.i+b.i)
    __radd__=__add__
    def __neg__(self): return C(-self.r,-self.i)
    def __sub__(self,b): return self+-b
    def __rsub__(self,b): return -self+b
    def __mul__(self,b):
        if not isinstance(b,C):return C(self.r*b,self.i*b)
        return C(self.r*b.r-self.i*b.i,self.r*b.i+self.i*b.r)
    __rmul__=__mul__
    def __truediv__(self,b):return C(self.r/b,self.i/b)
    def __eq__(self,b):return isinstance(b,C) and self.r==b.r and self.i==b.i
    def upper(self):return abs(self.r)+abs(self.i)
    def lower(self):return max(abs(self.r),abs(self.i))
    def wire(self):return [str(self.r),str(self.i)]
    def floating(self):return complex(float(self.r),float(self.i))


COEFS={(6,0):16,(5,0):-47,(4,1):-60,(4,0):-218,(3,2):-8,
       (3,1):-53,(3,0):-112,(2,2):68,(2,1):16,(2,0):84,
       (1,3):15,(1,2):26,(1,1):-84,(0,4):1,(0,3):-3,(0,2):21}
CERT=Path(__file__).resolve().parents[1]/'05-knowledge/results/planar_jc48_sep06_boundary_braid_certificate.json.gz'
GATES=0


def check(value,name):
    global GATES
    GATES+=1
    if not value:raise RuntimeError(name)


def powers(z,n):
    out=[C(1)]
    for k in range(n):out.append(out[-1]*z)
    return out


def base_taylor(u):
    up=powers(u,6);a=[[C() for k in range(5)] for j in range(7)]
    for (i,l),c in COEFS.items():
        for j in range(i+1):a[j][l]=a[j][l]+up[i-j]*(c*comb(i,j))
    return a


def taylor(a,z):
    zp=powers(z,4)
    return [[sum((a[j][l]*zp[l-k]*comb(l,k) for l in range(k,5)),C())
             for k in range(5)] for j in range(7)]


def rouche(b,r,du):
    rp=[r**k for k in range(5)];dp=[du**j for j in range(7)]
    lower=b[0][1].lower()*r
    upper=b[0][0].upper()+sum(b[0][k].upper()*rp[k] for k in range(2,5))
    upper+=sum(b[j][k].upper()*dp[j]*rp[k] for j in range(1,7) for k in range(5))
    return lower>upper


def radii(zs):
    return [min((zs[i]-zs[j]).lower() for j in range(4) if i!=j)/32 for i in range(4)]


def tube(u0,z0,u1,z1):
    rr=radii(z0);small=radii(z1)
    if any((z1[i]-z0[i]).upper()>=rr[i]/2 for i in range(4)):return False
    a=base_taylor(u0);end=base_taylor(u1);du=(u1-u0).upper()
    for i in range(4):
        if not rouche(taylor(a,z0[i]),rr[i],du):return False
        epsilon=min(rr[i],small[i])/64
        if not rouche(taylor(end,z1[i]),epsilon,Q(0)):return False
    return True


def braid_word(points):
    # Rational complex multiplier 1+i/4 sets a generic, oriented projection.
    project=lambda z:z*C(1,Q(1,4))
    first=[project(z) for z in points[0]]
    order=sorted(range(4),key=lambda j:first[j].r)
    check(order==list(range(4)),'basepoint root order')
    word=[]
    for raw0,raw1 in zip(points,points[1:]):
        z0,z1=[project(z) for z in raw0],[project(z) for z in raw1]
        events=[]
        for i in range(4):
            for j in range(i):
                a=(z0[i]-z0[j]).r;b=(z1[i]-z1[j]).r
                check(a!=0 and b!=0,'generic endpoint projection')
                if a*b<0:events.append((a/(a-b),i,j))
        check(len({e[0] for e in events})==len(events),'distinct crossing times')
        for h,i,j in sorted(events):
            ii,jj=order.index(i),order.index(j)
            check(abs(ii-jj)==1,'adjacent crossing')
            k=min(ii,jj);left,right=order[k:k+2]
            imag=((z0[left]-z0[right])*(1-h)+(z1[left]-z1[right])*h).i
            check(imag!=0,'nonintersecting polygonal strands')
            letter=(k+1)*(1 if imag<0 else -1)
            if word and word[-1]==-letter:word.pop()
            else:word.append(letter)
            order[k:k+2]=[right,left]
    return word


def produce():
    import numpy as np
    pp=list(permutations(range(4)))
    def quant(z):return C(Q(round(z.real*2**38),2**38),Q(round(z.imag*2**38),2**38))
    def roots(u):
        ts=np.roots([1,1,1,0,-u.floating()])
        return [quant(2*t**6+3*t**5+2*t**3+2*t**2) for t in ts]
    base=C(2,1)
    zbase=sorted(roots(base),key=lambda z:(z*C(1,Q(1,4))).r)
    centers=[C(),C(Q(-22528,2**20),Q(225902,2**20)),
             C(Q(-22528,2**20),Q(-225902,2**20)),
             C(Q(-612673,2**20),Q(2796037,2**20))]
    out=[]
    for index,center in enumerate(centers):
        r=Q(1,16) if index in (0,3) else Q(1,32)
        vertices=[base,center+C(r),center+C(0,r),center-C(r),center-C(0,r),center+C(r),base]
        us=[base];zs=[zbase];failures=0
        for start,end in zip(vertices,vertices[1:]):
            h=Q(1,16);done=Q(0)
            while done<1:
                h=min(h,1-done);u1=start+(end-start)*(done+h)
                raw=roots(u1);old=zs[-1]
                match=min(pp,key=lambda perm:sum(abs(raw[perm[k]].floating()-old[k].floating())**2 for k in range(4)))
                new=[raw[match[k]] for k in range(4)]
                if tube(us[-1],old,u1,new):
                    us.append(u1);zs.append(new);done+=h;h=min(2*h,Q(1,8))
                else:
                    failures+=1;h/=2
                    if h<Q(1,2**28):raise RuntimeError(('step underflow',index,done))
        word=braid_word(zs)
        out.append(dict(u=[u.wire() for u in us],z=[[z.wire() for z in row] for row in zs],word=word))
        print('produced loop',index,'steps',len(us)-1,'halvings',failures,'word',word,flush=True)
        raw=(json.dumps(dict(status='RATIONAL PATH WITNESS; verify before use',loops=out),sort_keys=True,separators=(',',':'))+'\n').encode()
        CERT.write_bytes(gzip.compress(raw,mtime=0))


def mul(a,b):return tuple(a[b[i]] for i in range(6))
def inv(a):return tuple(a.index(i) for i in range(6))


def group_census(words):
    cycles=[]
    for a,b,c in combinations(range(6),3):
        for q in ((a,b,c),(a,c,b)):
            p=list(range(6))
            for j in range(3):p[q[j]]=q[(j+1)%3]
            cycles.append(tuple(p))
    index={p:i for i,p in enumerate(cycles)}
    conjugate=[[index[mul(mul(a,b),inv(a))] for b in cycles] for a in cycles]
    inverse=[index[inv(a)] for a in cycles]
    rows=[(0,)+r for r in product(range(40),repeat=3)];counts=[]
    for word in words:
        survivors=[]
        for original in rows:
            row=list(original)
            for letter in word:
                k=abs(letter)-1;a,b=row[k:k+2]
                row[k:k+2]=[conjugate[a][b],a] if letter>0 else [b,conjugate[inverse[b]][a]]
            if tuple(row)==original:survivors.append(original)
        rows=survivors;counts.append(len(rows))
    transitive=[]
    for row in rows:
        orbit={0}
        while True:
            grow=orbit|{cycles[j][i] for j in row for i in orbit}
            if grow==orbit:break
            orbit=grow
        if len(orbit)==6:transitive.append(row)
    return counts,rows,transitive


def verify():
    compressed=CERT.read_bytes();raw=gzip.decompress(compressed);data=json.loads(raw);words=[];sizes=[]
    check(len(data['loops'])==4,'four declared actual loops')
    reference=None
    for number,loop in enumerate(data['loops']):
        us=[C(*p) for p in loop['u']]
        zs=[[C(*p) for p in row] for row in loop['z']]
        check(len(us)==len(zs),'path lengths agree')
        check(us[0]==us[-1]==C(2,1),'closed base loop')
        if reference is None:reference=zs[0]
        check(zs[0]==reference,'same labelled base fibre')
        check(sorted(z.wire() for z in zs[0])==sorted(z.wire() for z in zs[-1]),'closed unordered polygon fibre')
        start=base_taylor(us[0]);rr=radii(zs[0])
        for i in range(4):check(rouche(taylor(start,zs[0][i]),rr[i]/64,Q(0)),'initial root isolation')
        for k in range(len(us)-1):check(tube(us[k],zs[k],us[k+1],zs[k+1]),'exact disjoint Rouche tube')
        word=braid_word(zs);check(word==loop['word'],'exact polygon braid word')
        words.append(word);sizes.append(len(us)-1)
    counts,rows,transitive=group_census(words)
    check(not transitive,'no transitive six-sheet three-cycle tuples')
    print('RATIONAL ROUCHE CERTIFICATE PASS; analytic transport proof required separately')
    print('segment counts='+str(sizes))
    print('braid words='+str(words))
    print('40^3 fixed-first-cycle tuple counts='+str(counts))
    print('surviving tuples='+str(rows))
    print('transitive tuples='+str(transitive))
    print('certificate SHA256='+hashlib.sha256(raw).hexdigest())
    print('compressed certificate SHA256='+hashlib.sha256(compressed).hexdigest())
    print('PASS always-active gates='+str(GATES))


if __name__=='__main__':
    produce() if '--produce' in sys.argv else verify()
