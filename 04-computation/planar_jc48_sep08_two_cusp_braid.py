#!/usr/bin/env python3
"""Exact two-cusp rational configuration witnesses; proof in the companion note.

Default verifies four rational witnesses; --produce NAME proposes and certifies
one witness. --scout STEPS is heuristic and retains all six critical values.
"""
from fractions import Fraction as Q
from math import comb
from itertools import permutations
from pathlib import Path
import json,gzip,hashlib,sys

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


COEFS={}
CERT=None
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



def free_inverse(w):return tuple(-i for i in w[::-1])
def free_mul(*words):
    out=[]
    for w in words:
        for i in w:
            if out and out[-1]==-i:out.pop()
            else:out.append(i)
    return tuple(out)
def free_action(word,row=None):
    row=list(row or [(1,),(2,),(3,),(4,)])
    for i in word:
        j=abs(i)-1;a,b=row[j:j+2]
        row[j:j+2]=[free_mul(a,b,free_inverse(a)),a] if i>0 else[b,free_mul(free_inverse(b),a,b)]
    return row



COEFS={(6,0):Q(1),(5,0):Q(3,2),(4,1):Q(6),(4,0):Q(91,48),
 (3,2):Q(-2),(3,1):Q(145,18),(3,0):Q(925,162),(2,2):Q(15,2),
 (2,1):Q(-1,3),(2,0):Q(337,54),(1,3):Q(-6),(1,2):Q(38,3),
 (1,1):Q(-2),(1,0):Q(17,9),(0,4):Q(1),(0,3):Q(-4),(0,2):Q(34,9)}
PATHS={
 'smooth':(C(0),[3,1,-3]),
 'cusps':(C(-1),[3,2,3,1,3,1,1,3,-2,-3]),
 'node3':(C(3),[3,3]),
 'node_left':(C(Q(-401031,2**20)),[3,2,2,-3]),
 'node_upper':(C(Q(-1372351,2**20),Q(825221,2**20)),[3,2,3,1,2,2,-1,-3,-2,-3]),
 'node_lower':(C(Q(-1372351,2**20),Q(-825221,2**20)),[3,-2,-3,-1,2,2,1,3,2,-3])}
SELECTED=['smooth','cusps','node3','node_left']
BASE=C(4,3)

def vertices(name):
 c=PATHS[name][0];r=Q(1,32)
 return [BASE,c+C(r),c+C(0,r),c-C(r),c-C(0,r),c+C(r),BASE]
def certpath(name):
 return Path(__file__).resolve().parents[1]/('05-knowledge/results/planar_jc48_sep08_two_cusp_braid_'+name+'_certificate.json.gz')
def proposals():
 import numpy as np
 def roots(u):
  ts=np.roots([1,0,-2,0,-u.floating()])
  zs=[t**6-1.5*t**4+t**3/3-t for t in ts]
  return [C(Q(round(z.real*2**38),2**38),Q(round(z.imag*2**38),2**38)) for z in zs]
 return roots

def produce(name):
 if name not in SELECTED:raise RuntimeError('undeclared production path')
 roots=proposals();pp=list(permutations(range(4)))
 zbase=sorted(roots(BASE),key=lambda z:(z*C(1,Q(1,4))).r)
 us=[BASE];zs=[zbase];failures=0
 vv=vertices(name)
 for edge,(start,end) in enumerate(zip(vv,vv[1:])):
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
    if h<Q(1,2**30):raise RuntimeError(('step underflow',name,edge,done))
  print('EXACT EDGE',name,edge,'total segments',len(us)-1,flush=True)
 w=braid_word(zs);check(w==PATHS[name][1],'produced exact word')
 data=dict(status='RATIONAL WITNESS; verify before use',name=name,
   u=[u.wire() for u in us],z=[[z.wire() for z in row] for row in zs],word=w)
 raw=(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n').encode()
 certpath(name).write_bytes(gzip.compress(raw,mtime=0))
 print('EXACT PRODUCED',name,'segments',len(us)-1,'halvings',failures,'word',w,flush=True)

def scout(steps):
 roots=proposals();pp=list(permutations(range(4)))
 zbase=sorted(roots(BASE),key=lambda z:(z*C(1,Q(1,4))).r)
 for name in PATHS:
  vv=vertices(name);zs=[zbase]
  for start,end in zip(vv,vv[1:]):
   for j in range(1,steps+1):
    raw=roots(start+(end-start)*Q(j,steps));old=zs[-1]
    match=min(pp,key=lambda perm:sum(abs(raw[perm[k]].floating()-old[k].floating())**2 for k in range(4)))
    zs.append([raw[match[k]] for k in range(4)])
  print('HEURISTIC ONLY',steps,name,braid_word(zs),flush=True)

def geometry():
 import sympy as S
 t,u,v=S.symbols('t u v');rat=lambda c:S.Rational(c.numerator,c.denominator)
 U=t**4-2*t*t;V=t**6-S.Rational(3,2)*t**4+t**3/3-t
 F=sum(rat(c)*u**i*v**j for(i,j),c in COEFS.items())
 check(S.expand(F-S.resultant(U-u,V-v,t))==0,'actual resultant coefficients')
 check(S.expand(F.subs({u:U,v:V}))==0,'actual substitution')
 check(S.Poly(F,v).LC()==1 and S.degree(F,v)==4,'monic quartet')
 H=324*u**3+972*u*u+1080*u+289
 disc=-S.Rational(256,531441)*u*(u-3)**2*(u+1)**6*H**2
 check(S.factor(S.discriminant(F,v)-disc)==0,'all six critical values retained')
 check(S.discriminant(H,u)==-59592250800,'three cubic values distinct')
 check(S.resultant(H,u*(u-3)*(u+1),u)==868900175,'six values are distinct')
 check(S.factor(F.subs(u,-1)-(6*v-1)**2*(6*v+7)**2/1296)==0,'co-projected cusp targets')
 check(S.factor(F.subs(u,0)-v*v*(9*v*v-36*v+34)/9)==0,'single smooth fold')
 return str(S.factor(disc))

def algebra():
 w=PATHS['smooth'][1]
 check(free_action(w)==free_action([1]),'smooth word is first half twist')
 check(PATHS['node3'][1]==[3,3],'real node square')
 check(PATHS['node_left'][1]==[3,2,2,-3],'other real node conjugate square')
 check(free_action([3])==[(1,),(2,),(3,4,-3),(3,)],'other node prefix')
 check(free_action([3,2])==[(1,),(2,3,4,-3,-2),(2,),(3,)],'common cusp prefix')
 standard=[3,2]+[1,1,1,3,3,3]+[-2,-3]
 check(free_action(PATHS['cusps'][1])==free_action(standard),'two disjoint cusp cubes')
 for name in SELECTED:
  check(free_action(PATHS[name][1],[(1,)]*4)==[(1,)]*4,'cyclic converse control')
 return 'analytic necessity: a=b=d=c; no finite-group census'

def on_edges(us,name):
 vv=vertices(name);edge=0;last=Q(0)
 for u in us[1:]:
  start,end=vv[edge:edge+2];delta=end-start
  h=(u.r-start.r)/delta.r if delta.r else(u.i-start.i)/delta.i
  check(last<h<=1 and u==start+delta*h,'ordered rational edge subdivision')
  last=h
  if h==1:edge+=1;last=Q(0)
 check(edge==6,'six edges completed')

def verify(name,reference=None):
 gate0=GATES;compressed=certpath(name).read_bytes();raw=gzip.decompress(compressed)
 data=json.loads(raw);check(data['name']==name,'declared witness name')
 us=[C(*p) for p in data['u']];zs=[[C(*p) for p in row] for row in data['z']]
 check(len(us)==len(zs),'point counts agree')
 check(us[0]==us[-1]==BASE,'closed base loop')
 on_edges(us,name)
 if reference is not None:check(zs[0]==reference,'common labelled base fibre')
 check(sorted(z.wire() for z in zs[0])==sorted(z.wire() for z in zs[-1]),'closed unordered endpoint fibre')
 aa=base_taylor(BASE);rr=radii(zs[0])
 for i in range(4):check(rouche(taylor(aa,zs[0][i]),rr[i]/64,Q(0)),'initial root isolation')
 for k in range(len(us)-1):check(tube(us[k],zs[k],us[k+1],zs[k+1]),'uniform rational root tube')
 w=braid_word(zs);check(w==data['word']==PATHS[name][1],'actual rational crossing word')
 return zs[0],dict(name=name,segments=len(us)-1,word=w,gates=GATES-gate0,
  compressed_bytes=len(compressed),raw_bytes=len(raw),compressed_sha256=hashlib.sha256(compressed).hexdigest(),raw_sha256=hashlib.sha256(raw).hexdigest())

def main():
 if '--produce' in sys.argv:produce(sys.argv[-1]);return
 if '--scout' in sys.argv:scout(int(sys.argv[-1]));return
 disc=geometry();conclusion=algebra();reference=None;rows=[]
 for name in SELECTED:
  reference,row=verify(name,reference);rows.append(row)
 print('FINITE-EXACT TWO-CUSP RATIONAL BRAID PASS; analytic consumer proved separately')
 print(json.dumps(dict(discriminant=disc,paths=rows,conclusion=conclusion),sort_keys=True,indent=2))
 print('PASS always-active gates='+str(GATES))
if __name__=='__main__':main()
