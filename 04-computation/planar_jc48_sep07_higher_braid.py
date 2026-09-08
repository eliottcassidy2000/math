#!/usr/bin/env python3
"""Exact higher-cusp geometry and two-loop rational braid certificates.

Default: exact verification of all three literal curves and frozen witnesses.
--produce CASE: numerical proposals accepted only through exact rational tubes.
--scout CASE STEPS: heuristic words only, with no mathematical promotion.
No permutation census is used: the meridian reduction is checked in a free group.
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



CASES={
 'finite7_infinity9':{
  'v':[0,0,Q(-3,2),Q(-3,2),0,3,2], 'finite':7,'infinity':9,'nodes':3,
  'coefs':{(6,0):16,(5,0):-75,(4,1):-60,(4,0):Q(159,2),
   (3,2):-8,(3,1):Q(265,2),(3,0):Q(63,16),(2,2):75,(2,1):Q(93,8),
   (2,0):Q(63,16),(1,3):15,(1,2):Q(27,4),(1,1):Q(21,4),
   (0,4):1,(0,3):Q(1,2),(0,2):Q(7,4)},
  'words':[[3,2,3,2,1,-2,-3,-2,-3],[-3,-3,2,3,3]],
 },
 'finite7_infinity7':{
  'v':[0,0,0,0,0,0,2], 'finite':7,'infinity':7,'nodes':4,
  'coefs':{(6,0):16,(4,1):-72,(3,2):-8,(3,1):-16,(2,2):48,
   (0,4):1,(0,3):-4,(0,2):4},
  'words':[[2,-1,-2,1,2,-1,2,1,-2],[-1,-2,-3,-1,-1,2,1,1,3,2,1]],
 },
 'finite9_infinity7':{
  'v':[0,0,Q(-6,7),Q(-6,7),0,Q(12,7),2], 'finite':9,'infinity':7,'nodes':3,
  'coefs':{(6,0):16,(5,0):Q(-83568,2401),(4,1):Q(-2328,49),
   (4,0):Q(42384,2401),(3,2):-8,(3,1):Q(320,7),(3,0):Q(-3312,2401),
   (2,2):Q(1812,49),(2,1):Q(-120,49),(2,0):Q(1872,2401),
   (1,3):Q(60,7),(1,2):Q(-108,49),(1,1):Q(624,343),
   (0,4):1,(0,3):Q(-10,7),(0,2):Q(52,49)},
  'words':[[2,3,2,2,-1,2,1,-2,-2,-3,-2],
           [2,3,-1,-2,-1,-3,-1,2,1,3,1,2,1,-3,-2]],
 }
}


def configure(name):
    global COEFS,CERT
    if name not in CASES:raise RuntimeError('undeclared curve')
    COEFS=CASES[name]['coefs']
    CERT=Path(__file__).resolve().parents[1]/('05-knowledge/results/planar_jc48_sep07_higher_braid_'+name+'_certificate.json.gz')


def loop_centers():
    return [C(Q(-22528,2**20),Q(225902,2**20)),
            C(Q(-22528,2**20),Q(-225902,2**20))]


def vertices(center):
    base=C(2,1);r=Q(1,32)
    return [base,center+C(r),center+C(0,r),center-C(r),center-C(0,r),center+C(r),base]


def proposals(name):
    import numpy as np
    def roots(u):
        ts=np.roots([1,1,1,0,-u.floating()])
        return [sum(float(c)*t**k for k,c in enumerate(CASES[name]['v'])) for t in ts]
    def quant(z):return C(Q(round(z.real*2**38),2**38),Q(round(z.imag*2**38),2**38))
    return lambda u:[quant(z) for z in roots(u)]


def produce(name):
    configure(name);roots=proposals(name);pp=list(permutations(range(4)))
    base=C(2,1);zbase=sorted(roots(base),key=lambda z:(z*C(1,Q(1,4))).r)
    out=[]
    for index,center in enumerate(loop_centers()):
        us=[base];zs=[zbase];failures=0
        for start,end in zip(vertices(center),vertices(center)[1:]):
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
                    if h<Q(1,2**30):raise RuntimeError(('step underflow',name,index,done))
        word=braid_word(zs)
        out.append(dict(u=[u.wire() for u in us],z=[[z.wire() for z in row] for row in zs],word=word))
        print('EXACT PRODUCED',name,'loop',index,'steps',len(us)-1,'halvings',failures,'word',word,flush=True)
        raw=(json.dumps(dict(status='RATIONAL PATH WITNESS; verify before use',case=name,loops=out),sort_keys=True,separators=(',',':'))+'\n').encode()
        CERT.write_bytes(gzip.compress(raw,mtime=0))


def scout(name,steps):
    configure(name);roots=proposals(name);pp=list(permutations(range(4)))
    base=C(2,1);zbase=sorted(roots(base),key=lambda z:(z*C(1,Q(1,4))).r)
    for index,center in enumerate(loop_centers()):
        zs=[zbase]
        for start,end in zip(vertices(center),vertices(center)[1:]):
            for j in range(1,steps+1):
                raw=roots(start+(end-start)*Q(j,steps));old=zs[-1]
                match=min(pp,key=lambda perm:sum(abs(raw[perm[k]].floating()-old[k].floating())**2 for k in range(4)))
                zs.append([raw[match[k]] for k in range(4)])
        print('HEURISTIC ONLY',name,index,steps,braid_word(zs),flush=True)


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


def algebra(name,words):
    check(words==CASES[name]['words'],'literal certified words match analytic reduction')
    x,y=(1,),(2,)
    if name=='finite7_infinity9':
        upper=free_action([3,2,3,2]);lower=free_action([-3,-3])
        check(upper[0]==(1,) and upper[1]==(2,3,4,3,-4,-3,-2),'first arbitrary-group prefix')
        check(lower[1]==(2,) and lower[2]==(-4,3,4),'second arbitrary-group prefix')
        b=free_mul(free_inverse(y),x,y)
        a=free_mul(b,x,y,x,free_inverse(y),free_inverse(x),free_inverse(b))
        row=[a,b,x,y]
    else:
        prefix=[2,-1,-2,1] if name=='finite7_infinity7' else [2,2,-1]
        upper=free_action(prefix)
        if name=='finite7_infinity7':
            check(upper[1]==(2,3,-2) and upper[2]==(-3,-2,1,2,3),'seventh upper equality')
        else:
            check(upper[1]==(2,3,-2,-3,-2,1,2,3,2,-3,-2) and upper[2]==(2,3,-2),'ninth upper equality')
            check(free_action([2,3])==[(1,),(2,3,-2),(2,4,-2),(2,)],'common positive-meridian basis')
            for w,short in zip(words,[[2,2,-1,2,1,-2,-2],[-1,-2,-1,-3,-1,2,1,3,1,2,1]]):
                check(w==[2,3]+short+[-3,-2],'literal common conjugator')
        lower=free_action([-1,-2,-3,-1,-1] if name=='finite7_infinity7' else [-1,-2,-1,-3,-1])
        check(lower[1]==(-3,-2,3,2,3) and lower[2]==(4,),'last generator is eliminated')
        a=free_mul(x,y,x,y,free_inverse(x),free_inverse(y),free_inverse(x))
        d=free_mul(free_inverse(y),free_inverse(x),y,x,y)
        row=[a,x,y,d]
        if name=='finite9_infinity7':row=free_action([-3,-2],row)
    for w in words:check(free_action(w,row)==row,'free two-generator converse control')
    return 'two positive meridians; no finite-group enumeration'


def geometry(name):
    import sympy as S
    s,t,p,q,u,v,A,B,z=S.symbols('s t p q u v A B z')
    rat=lambda c:S.Rational(Q(c).numerator,Q(c).denominator)
    U=t**4+t**3+t**2;V=sum(rat(c)*t**k for k,c in enumerate(CASES[name]['v']))
    F=sum(rat(c)*u**i*v**j for (i,j),c in COEFS.items())
    def eq(a,b,label):check(S.cancel(a-b)==0,label)
    eq(F,S.resultant(U-u,V-v,t),'literal actual resultant')
    eq(F.subs({u:U,v:V}),0,'literal actual substitution')
    check(S.Poly(F,v).LC()==1 and S.degree(F,v)==4,'monic actual quartet')
    check(S.gcd(S.diff(U,t),S.diff(V,t))==t,'only critical normalization point')
    check(S.gcd(U,V)==t*t,'cusp has no second preimage')
    W=V
    for j in range(1,(CASES[name]['finite']+1)//2):
        W=S.expand(W-S.expand(W).coeff(t,2*j)*U**j)
    first=min(k[0] for k,c in S.Poly(W,t).terms())
    check(first==CASES[name]['finite'],'first unremoved odd finite coefficient')
    expected={'finite7_infinity9':S.Rational(9,2),'finite7_infinity7':-6,'finite9_infinity7':-S.Rational(44,7)}[name]
    eq(S.expand(W).coeff(t,first),expected,'finite leading coefficient')
    N=S.cancel((U.subs(t,s)-U)/(s-t));M=S.cancel((V.subs(t,s)-V)/(s-t))
    qp=p*(p*p+p+1)/(2*p+1)
    eq(S.rem(N.subs(t,p-s),s*s-p*s+q,s),p*(p*p+p+1)-(2*p+1)*q,'complete first pair equation')
    eq((p*(p*p+p+1)).subs(p,-S.Rational(1,2)),-S.Rational(3,8),'pair denominator is a unit')
    H={'finite7_infinity9':2*p**3+6*p*p+10*p+9,
       'finite7_infinity7':(p-1)*(p+1)*(p*p+2*p+3),
       'finite9_infinity7':7*p**3+20*p*p+32*p+22}[name]
    expected_pair={'finite7_infinity9':-p**3*H/(2*(2*p+1)),
       'finite7_infinity7':-2*p**3*H/(2*p+1)**2,
       'finite9_infinity7':-2*p**4*H/(7*(2*p+1)**2)}[name]
    eq(S.rem(M.subs(t,p-s),s*s-p*s+q,s).subs(q,qp),expected_pair,'complete residual pair polynomial')
    eq(p*p-4*qp,-p*(2*p*p+3*p+4)/(2*p+1),'pair discriminant')
    check(S.discriminant(H,p)!=0,'all node pairs distinct')
    check(S.resultant(H,p*(2*p+1)*(2*p*p+3*p+4),p)!=0,'pairs are off diagonal')
    check(S.degree(H,p)==CASES[name]['nodes'],'declared number of node pairs')
    T=S.diff(U,t).subs(t,s)*S.diff(V,t)-S.diff(V,t).subs(t,s)*S.diff(U,t)
    gb=[s-2*t**5+t**3+t*t+t,t**6] if CASES[name]['finite']==7 else [s-11*t**7-6*t**6-2*t**5+t**3+t*t+t,t**8]
    check(S.groebner([N,M,T],s,t,domain=S.QQ)==S.groebner(gb,s,t,domain=S.QQ),'every off-diagonal pair transverse')
    R=S.rem(V-B,U-A,t)
    check(S.degree(R,t)==3 and not S.Poly(R,t).LC().free_symbols,'triple remainder has constant nonzero leading term')
    coeffs=[S.together(c).as_numer_denom()[0] for c in S.Poly(S.rem(U-A,R,t),t).all_coeffs()]
    check(S.groebner(coeffs,A,B,domain=S.QQ)==S.groebner([1],A,B,domain=S.QQ),'no triple image')
    X=S.cancel(U.subs(t,1/z)/V.subs(t,1/z));Z=S.cancel(1/V.subs(t,1/z))
    eq(S.limit(X/z**2,z,0),S.Rational(1,2),'infinity first coordinate order')
    eq(S.limit(Z/z**6,z,0),S.Rational(1,2),'line contact six')
    if CASES[name]['infinity']==9:
        eq(S.limit((Z-4*X**3)/z**7,z,0),0,'seventh infinity coefficient vanishes')
        eq(S.limit((Z-4*X**3)/X**4,z,0),-30,'even infinity coefficient')
        eq(S.limit((Z-4*X**3+30*X**4)/z**9,z,0),S.Rational(7,16),'first odd infinity ninth coefficient')
    else:
        expected=-S.Rational(3,2) if name=='finite7_infinity7' else -S.Rational(9,14)
        eq(S.limit((Z-4*X**3)/z**7,z,0),expected,'first odd infinity seventh coefficient')
    check((CASES[name]['finite']-1)//2+(CASES[name]['infinity']-1)//2+CASES[name]['nodes']==10,'complete sextic genus ledger')
    disc={'finite7_infinity9':-u**7*(256*u*u+11*u+12)*(3136*u**3+6032*u*u+4452*u+1323)**2/4096,
      'finite7_infinity7':-4096*u**9*(u+2)**2*(9*u*u-4*u+6)**2*(256*u*u+11*u+12),
      'finite9_infinity7':-S.Rational(4096,7**12)*u**9*(256*u*u+11*u+12)*(194481*u**3+87242*u*u+62998*u+11154)**2}[name]
    eq(S.discriminant(F,v),disc,'complete actual vertical discriminant')
    if name=='finite7_infinity7':
        eq(F.subs(u,0),v*v*(v-2)**2,'co-projected cusp and node retained')
        eq(S.rem(U,t*t+t+1,t),0,'node projection agrees with cusp projection')
        eq(S.rem(V,t*t+t+1,t),2,'node target differs from cusp target')
    return dict(finite=CASES[name]['finite'],infinity=CASES[name]['infinity'],nodes=CASES[name]['nodes'],pair=str(S.expand(H)),discriminant=str(S.factor(disc)))


def on_declared_edges(us,center):
    vv=vertices(center);edge=0;last=Q(0)
    for u in us[1:]:
        start,end=vv[edge:edge+2];delta=end-start
        h=(u.r-start.r)/delta.r if delta.r else (u.i-start.i)/delta.i
        check(last<h<=1 and u==start+delta*h,'exact ordered rational edge subdivision')
        last=h
        if h==1:edge+=1;last=Q(0)
    check(edge==6,'all six declared edges completed')


def verify_case(name):
    configure(name);start_gate=GATES;geo=geometry(name)
    compressed=CERT.read_bytes();raw=gzip.decompress(compressed);data=json.loads(raw)
    check(data['case']==name,'literal witness curve identity')
    check(len(data['loops'])==2,'two declared actual loops')
    words=[];sizes=[];reference=None
    for number,loop in enumerate(data['loops']):
        us=[C(*p) for p in loop['u']];zs=[[C(*p) for p in row] for row in loop['z']]
        check(len(us)==len(zs),'path lengths agree')
        check(us[0]==us[-1]==C(2,1),'closed common base loop')
        on_declared_edges(us,loop_centers()[number])
        if reference is None:reference=zs[0]
        check(zs[0]==reference,'same labelled base fibre')
        check(sorted(z.wire() for z in zs[0])==sorted(z.wire() for z in zs[-1]),'closed unordered polygon fibre')
        aa=base_taylor(us[0]);rr=radii(zs[0])
        for i in range(4):check(rouche(taylor(aa,zs[0][i]),rr[i]/64,Q(0)),'initial root isolation')
        for k in range(len(us)-1):check(tube(us[k],zs[k],us[k+1],zs[k+1]),'uniform disjoint rational Rouche tube')
        word=braid_word(zs);check(word==loop['word'],'exact rational braid word')
        words.append(word);sizes.append(len(us)-1)
    conclusion=algebra(name,words)
    return dict(case=name,geometry=geo,segments=sizes,words=words,conclusion=conclusion,
      gates=GATES-start_gate,compressed_bytes=len(compressed),raw_bytes=len(raw),
      compressed_sha256=hashlib.sha256(compressed).hexdigest(),raw_sha256=hashlib.sha256(raw).hexdigest())


def main():
    if '--produce' in sys.argv:
        produce(sys.argv[-1]);return
    if '--scout' in sys.argv:
        scout(sys.argv[-2],int(sys.argv[-1]));return
    selected=[sys.argv[-1]] if '--case' in sys.argv else list(CASES)
    reports=[verify_case(name) for name in selected]
    print('EXACT HIGHER-CUSP ROUCHE CERTIFICATE PASS; analytic transport proved separately')
    print(json.dumps(reports,sort_keys=True,indent=2))
    print('PASS always-active gates='+str(GATES))


if __name__=='__main__':main()
