#!/usr/bin/env python3
"""Exact mixed (5,3,3)-cusp witnesses; companion proof and audit required.

Default verifies six rational moving-centre witnesses; --produce NAME proposes and certifies
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


def affine_distance(d0,delta):
    # Min over h in [0,1] of max(|Re(d0+h delta)|,|Im(d0+h delta)|).
    candidates={Q(0),Q(1)}
    for a,b in [(d0.r,delta.r),(d0.i,delta.i),
                (d0.r-d0.i,delta.r-delta.i),(d0.r+d0.i,delta.r+delta.i)]:
        if b:
            h=-a/b
            if 0<=h<=1:candidates.add(h)
    return min((d0+delta*h).lower() for h in candidates)


def moving_taylor(b,du,dz):
    # F(u0+h du,z0+h dz+w), coefficient of h^j w^k.
    up=powers(du,6);zp=powers(dz,4);out={}
    for j in range(7):
        for k in range(5):
            if b[j][k]==C():continue
            for ell in range(k+1):
                key=(j+k-ell,ell)
                out[key]=out.get(key,C())+b[j][k]*up[j]*zp[k-ell]*comb(k,ell)
    return out


def moving_margin(coeffs,r):
    lhs=coeffs.get((0,1),C()).lower()*r
    rhs=sum((v.upper()*r**k for (j,k),v in coeffs.items() if (j,k)!=(0,1)),Q(0))
    return lhs-rhs


def tube(u0,z0,u1,z1):
    distances={(i,j):affine_distance(z0[i]-z0[j],z1[i]-z0[i]-z1[j]+z0[j])
               for i in range(4) for j in range(i)}
    if not all(distances.values()):return False
    rr=[min(d for(i,j),d in distances.items() if k in(i,j))/16 for k in range(4)]
    if any(rr[i]+rr[j]>=d for(i,j),d in distances.items()):return False
    aa=base_taylor(u0)
    return all(moving_margin(moving_taylor(taylor(aa,z0[i]),u1-u0,z1[i]-z0[i]),rr[i])>0
               for i in range(4))


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



COEFS={(6, 0):Q(1,1),
 (5, 0):Q(259,2),
 (4, 1):Q(-46,3),
 (4, 0):Q(499099,48),
 (3, 2):Q(-2,1),
 (3, 1):Q(-92903,18),
 (3, 0):Q(-1426528,27),
 (2, 2):Q(44885,54),
 (2, 1):Q(587489,27),
 (2, 0):Q(672224,9),
 (1, 3):Q(-146,3),
 (1, 2):Q(-215668,81),
 (1, 1):Q(-658112,27),
 (1, 0):Q(-67072,3),
 (0, 4):Q(1,1),
 (0, 3):Q(61532,729),
 (0, 2):Q(159428,81),
 (0, 1):Q(33536,9)}
PATHS={
 'cusp_plus':(C(Q(13,3),Q(0,1)),[-1, -2, -3, 2, 2, 2, 3, 2, 1]),
 'cusp_five':(C(Q(-19,3),Q(0,1)),[3, 2, 1, 3, 3, 3, 3, 3, -1, -2, -3]),
 'cusp_two':(C(Q(8,3),Q(0,1)),[3, 2, 1, 1, 1, -2, -3]),
 'node0':(C(Q(3077421,1048576),Q(0,1)),[-1, 3, -2, 1, 2, 2, -1, 2, -3, 1]),
 'node1':(C(Q(3,1),Q(0,1)),[-1, 3, -2, 1, 1, 2, -3, 1]),
 'node2':(C(Q(904195,262144),Q(0,1)),[-1, 3, -2, 3, 3, 2, -3, 1])}
SELECTED=['cusp_plus', 'cusp_five', 'cusp_two', 'node0', 'node1', 'node2']
LOCAL={'cusp_plus': {'centres': [['63/2', '0'], ['2944373252607/68719476736', '0'], ['5672540810507/274877906944', '0']], 'radii': ['1/64', '1/64', '1/64'], 'slopes': [['21/2', '0'], ['23908014653/1073741824', '0'], ['2899421239/536870912', '0']], 'pair_positions': [1, 2], 'stem': [-1, -2, -3]}, 'cusp_five': {'centres': [['-65/2', '0'], ['-45028884891915/274877906944', '-16244794860909/274877906944'], ['-45028884891915/274877906944', '16244794860909/274877906944']], 'radii': ['1/64', '1/64', '1/64'], 'slopes': [['9/2', '0'], ['21295879509/1073741824', '1265416875/134217728'], ['21295879509/1073741824', '-1265416875/134217728']], 'pair_positions': [2, 3], 'stem': [3, 2, 1]}, 'cusp_two': {'centres': [['8', '0'], ['3226007533537/274877906944', '0'], ['4847635461641/274877906944', '0']], 'radii': ['1/32', '1/64', '1/64'], 'slopes': [['18', '0'], ['708567703/134217728', '0'], ['3966094073/536870912', '0']], 'pair_positions': [0, 1], 'stem': [3, 2, -1]}, 'node0': {'centres': [['3615921682829/274877906944', '0'], ['675004297611/34359738368', '0'], ['3427503766705/274877906944', '0']], 'radii': ['5/32', '1/64', '1/64'], 'slopes': [['13479737067/1073741824', '0'], ['2040999969/268435456', '0'], ['17131898293/1073741824', '0']], 'pair_positions': [1, 2], 'stem': [-1, 3, -2, 1]}, 'node1': {'centres': [['27/2', '0'], ['3972431054581/274877906944', '0'], ['2768342727319/137438953472', '0']], 'radii': ['7/64', '1/64', '1/64'], 'slopes': [['21/2', '0'], ['10741043761/536870912', '0'], ['8224769609/1073741824', '0']], 'pair_positions': [0, 1], 'stem': [-1, 3, -2]}, 'node2': {'centres': [['13018346298397/549755813888', '0'], ['4367822953393/274877906944', '0'], ['21696011375/1073741824', '0']], 'radii': ['15/64', '1/64', '1/64'], 'slopes': [['15630568517/1073741824', '0'], ['716284629/134217728', '0'], ['15264028763/1073741824', '0']], 'pair_positions': [2, 3], 'stem': [-1, 3, -2, -3]}}

BASE=C(4,3)

def vertices(name):
 c=PATHS[name][0];r=Q(1,64) if name in ('node0','node1') else Q(1,32)
 return [BASE,c+C(r),c+C(0,r),c-C(r),c-C(0,r),c+C(r),BASE]
def certpath(name):
 return Path(__file__).resolve().parents[1]/('05-knowledge/results/planar_jc48_sep08_mixed_cusp_braid_'+name+'_certificate.json.gz')
def proposals():
 import numpy as np
 def roots(u):
  ts=np.roots([1,-8/3,-2,8,-u.floating()])
  zs=[t**6-1.5*t**4-16*t**3+48*t for t in ts]
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
 t,u,v,z,h=S.symbols('t u v z h');rat=lambda c:S.Rational(c.numerator,c.denominator)
 U=t**4-S.Rational(8,3)*t**3-2*t*t+8*t
 V=t**6-S.Rational(3,2)*t**4-16*t**3+48*t
 F=sum(rat(c)*u**i*v**j for(i,j),c in COEFS.items())
 check(S.expand(F-S.resultant(U-u,V-v,t))==0,'actual resultant coefficients')
 check(S.expand(F.subs({u:U,v:V}))==0,'actual substitution')
 check(S.Poly(F,v).LC()==1 and S.degree(F,v)==4,'monic quartet')
 check(S.Poly(F,u,v).total_degree()==6,'actual projective sextic')
 H=(u-3)*(2187*u*u-13962*u+22139)
 critical=(3*u-13)*(3*u-8)*(3*u+19)
 disc=-S.Rational(1048576,847288609443)*(3*u-13)**3*(3*u-8)**3*(3*u+19)**5*H**2
 check(S.factor(S.discriminant(F,v)-disc)==0,'all six critical values retained')
 check(S.factor(S.discriminant(U-u,t)+S.Rational(256,27)*critical)==0,'parameter discriminant')
 check(S.discriminant(H,u)!=0,'three node values distinct')
 check(S.resultant(H,critical,u)!=0,'six values are distinct')
 check(S.Poly(S.gcd(S.diff(U,t),S.diff(V,t)),t).monic().as_expr()==S.expand((t*t-1)*(t-2)),'exact three critical locus')
 for e in [-1,1,2]:
  check(S.Poly(S.gcd(U-U.subs(t,e),V-V.subs(t,e)),t).monic().as_expr()==S.expand((t-e)**2),'single cusp preimage')
  jet=(S.diff(U,t,2)*S.diff(V,t,3)-S.diff(V,t,2)*S.diff(U,t,3)).subs(t,e)
  check((jet==0)==(e==-1),'only prescribed fifth cusp has zero ordinary jet')
 uu=S.expand(U.subs(t,h-1)-U.subs(t,-1))
 ww=S.expand(V.subs(t,h-1)-V.subs(t,-1)-S.Rational(9,2)*uu-S.Rational(1,16)*uu**2)
 check(all(S.Poly(ww,h).nth(j)==0 for j in range(5)) and S.Poly(ww,h).nth(5)==4,'finite cusp exactly two-five')
 check(S.diff(U,t,3).subs(t,-1)!=0,'quadratic deck fixed-centre hostile excluded')
 def lead(expr,k,value,label):
  num,den=S.cancel(expr).as_numer_denom();poly=S.Poly(num,z)
  check(all(j[0]>=k for j,c in poly.terms()) and den.subs(z,0)!=0,label+' order')
  check(poly.coeff_monomial(z**k)==value*den.subs(z,0),label+' coefficient')
 X=S.cancel(U.subs(t,1/z)/V.subs(t,1/z));Z=S.cancel(1/V.subs(t,1/z))
 lead(X,2,1,'infinity transverse coordinate')
 lead(Z-X**3,7,8,'infinity odd type seven')
 check(2+1+1+3+3==10,'rational sextic genus budget')
 return str(S.factor(disc))

def algebra():
 def inverse(w):return [-i for i in w[::-1]]
 def conjugate(h,w):return free_mul(h,w,free_inverse(h))
 splits={'cusp_plus':([-1,-2,-3],[2,2,2],1),'cusp_five':([3,2,1],[3,3,3,3,3],2),
         'cusp_two':([3,2],[1,1,1],0),'node0':([-1,3,-2,1],[2,2],1),
         'node1':([-1,3,-2],[1,1],0),'node2':([-1,3,-2],[3,3],2)}
 generators=[(i,) for i in range(1,5)]
 x,y,z,w=generators;P0=[-1,3,-2]
 old=free_action(inverse(P0))
 expected_old=[free_mul(x,y,z,free_inverse(y),free_inverse(x)),x,w,conjugate(free_inverse(w),y)]
 check(old==expected_old,'exact old generators in positive new free basis')
 check(free_action(P0,old)==generators,'free change of basis has explicit inverse')
 check(free_action(inverse(P0),free_action(P0))==generators,'free inverse in both directions')
 h=free_mul(x,y)
 pairs=[(w,conjugate(free_inverse(w),y)),(x,w),(conjugate(h,z),conjugate(h,y)),(x,z),(x,y),(z,w)]
 conjugators=[w,(),free_inverse(h),(),(),()]
 h4_pairs=[(w,y),(x,w),(z,y),(x,z),(x,y),(z,w)]
 for (name,(prefix,core,k)),pair,g,target in zip(splits.items(),pairs,conjugators,h4_pairs):
  check(PATHS[name][1]==prefix+core+inverse(prefix),'literal conjugated word decomposition')
  check(tuple(free_action(prefix,old)[k:k+2])==pair,'exact original core pair in new free basis')
  check(tuple(conjugate(g,j) for j in pair)==target,'exact simultaneous pair transport to H4')
 # H^odd fixes a pair iff its second coordinate is fixed; product is invariant.
 a,b=(1,),(2,)
 for odd in (3,5):
  action=free_action([1]*odd,[a,b,(3,),(4,)])[:2]
  g=free_mul(*([a,b]*((odd-1)//2)))
  check(action[1]==conjugate(g,a),'odd Hurwitz second coordinate')
  check(free_mul(*action)==free_mul(a,b),'Hurwitz preserves pair product')
  check(free_mul(g,a)==tuple([1,2]*(odd//2)+[1]),'odd left alternating word')
  check(free_mul(b,g)==tuple([2,1]*(odd//2)+[2]),'odd right alternating word')
 # Actual local basis corrections: only the third cusp and final node acquire H^-.
 for name,(prefix,core,k) in splits.items():
  extra=[-(k+1)] if name in ('cusp_two','node2') else []
  check(LOCAL[name]['stem']==prefix+extra,'actual stem differs only inside the declared pair')
  check(LOCAL[name]['pair_positions']==[k,k+1],'actual adjacent local pair positions')
 return 'six fixed-word relations are exactly Artin H4 in a positive free-meridian basis; actual complement only a quotient; no Keller exclusion inferred'


def method_controls():
 import sympy as S
 h,w,u,v=S.symbols('h w u v');to_s=lambda z:S.Rational(z.r.numerator,z.r.denominator)+S.I*S.Rational(z.i.numerator,z.i.denominator)
 F=sum(S.Rational(c.numerator,c.denominator)*u**i*v**j for(i,j),c in COEFS.items())
 u0=C(Q(2,3),Q(1,5));v0=C(Q(-3,7),Q(4,9));du=C(Q(5,11),Q(-2,13));dv=C(Q(-7,17),Q(3,19))
 actual=moving_taylor(taylor(base_taylor(u0),v0),du,dv)
 literal=S.Poly(S.expand(F.subs({u:to_s(u0)+h*to_s(du),v:to_s(v0)+h*to_s(dv)+w})),h,w)
 check(all(S.expand(to_s(actual.get(key,C()))-literal.coeff_monomial(h**key[0]*w**key[1]))==0
           for key in set(actual)|set(literal.monoms())),'independent symbolic affine substitution')
 check(affine_distance(C(1),C(-2))==0,'swapped endpoint collision hostile')
 check(affine_distance(C(1),C(-1,1))==Q(1,2),'interior minimum retained')
 check(affine_distance(C(2,-1),C(0,2))==2,'flat interior maximum control')
 # Vanishing linear residual is insufficient: w-h^2+h^3 moves its root.
 check(moving_margin({(0,1):C(1),(2,0):C(-1),(3,0):C(1)},Q(1,16))<0,'nonlinear curvature hostile')
 # Unbounded common translation cancels exactly for F=(v-u)...(v-u-3).
 F0=S.prod(v-u-j for j in range(4));delta=C(10**6,10**6)
 for root in range(4):
  lit=S.Poly(S.expand(F0.subs(v,v+root)),u,v)
  b=[[C(lit.coeff_monomial(u**j*v**k)) for k in range(5)] for j in range(7)]
  moved=moving_taylor(b,delta,delta)
  check(all(value==C() for(j,k),value in moved.items() if j>0),'exact cancellation of unbounded common translation')
  check(moving_margin(moved,Q(1,16))>0,'long translation remains certifiable')


def local_cluster(name,us,zs):
 # A complete complex-parameter disk, not only the diamond boundary.
 data=LOCAL[name];c=PATHS[name][0];ur=Q(1,64) if name in ('node0','node1') else Q(1,32)
 centres=[C(*p) for p in data['centres']];slopes=[C(*p) for p in data['slopes']]
 rr=[Q(r) for r in data['radii']];margins=[]
 for i in range(3):
  coeffs=moving_taylor(taylor(base_taylor(c),centres[i]),C(1),slopes[i]);degree=2 if i==0 else 1
  lhs=coeffs.get((0,degree),C()).lower()*rr[i]**degree
  rhs=sum((v.upper()*ur**j*rr[i]**k for(j,k),v in coeffs.items() if(j,k)!=(0,degree)),Q(0))
  check(lhs>rhs,'whole complex-parameter disk cluster Rouche')
  margins.append(str(lhs-rhs))
 for i in range(3):
  for j in range(i):
   check((centres[i]-centres[j]).lower()-(slopes[i]-slopes[j]).upper()*ur>rr[i]+rr[j],
         'uniform holomorphic cluster separation')
 # Exactly one discriminant point, also inside the smaller disk contained in the diamond.
 ht=[2187*c*c*c-20523*c*c+64025*c-66417,6561*c*c-41046*c+64025,6561*c-20523,C(2187)]
 if name.startswith('cusp'):
  check(ht[0].lower()>sum(ht[j].upper()*ur**j for j in range(1,4)),'no node value in cusp parameter disk')
 else:
  for radius in [ur,ur/2]:
   check(ht[1].lower()*radius>ht[0].upper()+sum(ht[j].upper()*radius**j for j in range(2,4)),
         'unique node value in parameter disk and inside diamond')
 for value in [C(Q(13,3)),C(Q(-19,3)),C(Q(8,3))]:
  check(value==c or (value-c).lower()>ur,'no other cusp in parameter disk')
 endpoint=us.index(vertices(name)[1]);end=zs[endpoint];epsilon=[min(r/64,Q(1,2**20)) for r in radii(end)]
 aa=base_taylor(us[endpoint]);discs=[centres[i]+slopes[i]*ur for i in range(3)];assignment=[]
 for k in range(4):
  check(rouche(taylor(aa,end[k]),epsilon[k],Q(0)),'stem endpoint tiny root isolation')
  fits=[j for j in range(3) if(end[k]-discs[j]).upper()+epsilon[k]<rr[j]]
  check(len(fits)==1,'endpoint root belongs to one certified local cluster')
  assignment.append(fits[0])
 check([assignment.count(j) for j in range(3)]==[2,1,1],'actual two-one-one local root partition')
 project=lambda z:z*C(1,Q(1,4));order=sorted(range(4),key=lambda k:project(end[k]).r)
 for i,j in zip(order,order[1:]):
  check(project(end[j]-end[i]).r>Q(5,4)*(epsilon[i]+epsilon[j]),'actual endpoint projection order')
 pair=[order.index(k) for k in range(4) if assignment[k]==0]
 check(sorted(pair)==data['pair_positions'] and abs(pair[0]-pair[1])==1,'actual cluster occupies declared adjacent meridian positions')
 for k in range(4):
  if assignment[k]!=0:
   check(abs(project(end[k]-discs[0]).r)>Q(5,4)*(rr[0]+epsilon[k]),
         'other punctures outside projected local disk strip')
 stem=braid_word(zs[:endpoint+1]);check(stem==data['stem'],'exact actual incoming stem braid')
 return dict(stem_word=stem,pair_positions=sorted(pair),cluster_degrees=[2,1,1],cluster_radii=data['radii'])

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
 check(all(len(row)==4 for row in zs),'every stored fibre has four labelled centres')
 check(us[0]==us[-1]==BASE,'closed base loop')
 on_edges(us,name)
 if reference is not None:check(zs[0]==reference,'common labelled base fibre')
 check(sorted(z.wire() for z in zs[0])==sorted(z.wire() for z in zs[-1]),'closed unordered endpoint fibre')
 aa=base_taylor(BASE);rr=radii(zs[0])
 for i in range(4):check(rouche(taylor(aa,zs[0][i]),rr[i]/64,Q(0)),'initial root isolation')
 for k in range(len(us)-1):check(tube(us[k],zs[k],us[k+1],zs[k+1]),'uniform rational root tube')
 w=braid_word(zs);check(w==data['word']==PATHS[name][1],'actual rational crossing word')
 local=local_cluster(name,us,zs)
 return zs[0],dict(name=name,segments=len(us)-1,word=w,local=local,gates=GATES-gate0,
  compressed_bytes=len(compressed),raw_bytes=len(raw),compressed_sha256=hashlib.sha256(compressed).hexdigest(),raw_sha256=hashlib.sha256(raw).hexdigest())

def main():
 if '--produce' in sys.argv:produce(sys.argv[-1]);return
 if '--scout' in sys.argv:scout(int(sys.argv[-1]));return
 disc=geometry();conclusion=algebra();method_controls();reference=None;rows=[]
 for name in SELECTED:
  reference,row=verify(name,reference);rows.append(row)
 print('FINITE-EXACT MIXED-CUSP RATIONAL BRAID PASS; marked Artin H4 quotient; no Keller exclusion claimed')
 print(json.dumps(dict(discriminant=disc,paths=rows,conclusion=conclusion),sort_keys=True,indent=2))
 print('PASS always-active gates='+str(GATES))
if __name__=='__main__':main()
