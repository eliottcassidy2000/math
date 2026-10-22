"""Dyadic outward interval certificate for rank88 scalar extrema; no floats in decisions."""
from fractions import Fraction as Q
from pathlib import Path
import json,time
BITS=180; S=1<<BITS
class I:
 def __init__(self,a,b=None,raw=False):
  if raw:self.l,self.h=a,b;return
  a=Q(a);b=Q(b) if b is not None else a
  self.l=a.numerator*S//a.denominator;self.h=-((-b.numerator*S)//b.denominator)
 @classmethod
 def raw(cls,l,h):return cls(l,h,True)
 def __add__(a,b):
  if not isinstance(b,I):b=I(b)
  return I.raw(a.l+b.l,a.h+b.h)
 __radd__=__add__
 def __neg__(a):return I.raw(-a.h,-a.l)
 def __sub__(a,b):return a+-b if isinstance(b,I) else a+I(-b)
 def __rsub__(a,b):return I(b)+-a
 def __mul__(a,b):
  if not isinstance(b,I):b=I(b)
  z=[a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h]
  return I.raw(min(z)//S,-((-max(z))//S))
 __rmul__=__mul__
 def inv(a):
  if a.l<=0<=a.h:raise ZeroDivisionError
  return I.raw(S*S//a.h,-((-S*S)//a.l))
 def __truediv__(a,b):return a*(b if isinstance(b,I) else I(b)).inv()
 def __pow__(a,n):
  if n<0:return (a**(-n)).inv()
  z=I(1)
  while n:
   if n&1:z=z*a
   n//=2
   if n:a=a*a
  return z
 def bounds(a):return [float(Q(a.l,S)),float(Q(a.h,S))]
def atan(z,N):
 p=sum(((-1)**k*z**(2*k+1)/(2*k+1) for k in range(N)),Q(0));e=(-1)**N*z**(2*N+1)/(2*N+1)
 return I(min(p,p+e),max(p,p+e))
PI=16*atan(Q(1,5),90)-4*atan(Q(1,239),30)
def expminus(x):
 if not (0 <= x.l <= x.h <= 16*S):raise ValueError("expminus certified only on [0,16]")
 z=x/16;v=I(1);su=v
 for j in range(1,61):v=v*z/j;su=su+v
 # z <= 1; remainder < 2/61! < 2^-270, safely below 2^-250.
 su=su+I(0,Q(1,2**250))
 return su**-16
assert expminus(2*PI).h*500 < S
K=80;N=K+4
mul=lambda a,b:[sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N+1)]
def power(a,n):
 z=[1]+[0]*N
 for _ in range(n):z=mul(z,a)
 return z
sig=lambda n,p:sum(d**p for d in range(1,n+1) if n%d==0)
e4=[1]+[240*sig(n,3) for n in range(1,N+1)]
e6=[1]+[-504*sig(n,5) for n in range(1,N+1)]
e2=[1]+[-24*sig(n,1) for n in range(1,N+1)]
delta=[(a-b)//1728 for a,b in zip(power(e4,3),power(e6,2))]
f=[0]*(N+1)
for j,c in enumerate([1,-2640,1813680,-244992000]):
 v=mul(power(delta,j),power(e4,11-3*j))
 f=[a+c*b for a,b in zip(f,v)]
assert f[:5]==[1,0,0,0,168498000]
u=[n*f[n] for n in range(N+1)];v=[n*n*f[n] for n in range(N+1)];w=[n**3*f[n] for n in range(N+1)]
e2u=mul(e2,u);e2v=mul(e2,v)
A=[-w[n]+4*e2v[n] for n in range(4,N+1)]
B=[45*(v[n]-4*e2u[n]) for n in range(4,N+1)]
C=[45*u[n] for n in range(4,N+1)]
# Cauchy majorant |f_n| <= M R^-n on |q|=R.
R=Q(1,100);E4=1+240*R/(1-16*R);E6=1+504*R/(1-64*R);DD=(E4**3+E6**2)/1728
M=sum(c*DD**j*E4**(11-3*j) for j,c in enumerate([1,2640,1813680,244992000]))
# q <= 1/500 throughout t>=1. Polynomial and differentiated tails.
# n=K+5; ratio at most (1/5)*(1+1/n)^4 < 1/4.
n=K+5
assert Q(1,5)*(Q(n+1,n))**4<Q(1,4)
errP=M*R**-4*Q(1,5)**(K+1)*2000*n**3*Q(4,3)
errD=M*R**-4*Q(1,5)**(K+1)*3000*n**4*Q(4,3)
ERR=Q(1,10**20)
assert errP<ERR and errD<ERR

def poly(cs,q):
 z=I(0)
 for c in reversed(cs):z=z*q+c
 return z

# Coefficient arrays for successive derivatives along x, q=exp(-x).
DER=[(A,B,C)]
for order in range(3):
 aa,bb,cc=DER[-1]
 DER.append(([-k*c for k,c in enumerate(aa)], [2*a-k*b for k,(a,b) in enumerate(zip(aa,bb))], [b-k*c for k,(b,c) in enumerate(zip(bb,cc))]))
# Higher derivative tail bound. At m<=3, 20000*n^6 dominates coefficients.
assert Q(1,5)*(Q(n+1,n))**6<Q(1,4)
errHigher=M*R**-4*Q(1,5)**(K+1)*20000*n**6*Q(4,3)
assert errHigher<ERR

def rawstation(t,order=0):
 x=2*PI*t;q=expminus(x);aa,bb,cc=DER[order]
 return x*x*poly(aa,q)+x*poly(bb,q)+poly(cc,q)+I(-ERR,ERR)
def station(t,derivative=False):
 order=int(derivative)
 if t.l==t.h:return rawstation(t,order)
 mid=(t.l+t.h)//2;center=I.raw(mid,mid)
 return rawstation(center,order)+rawstation(t,order+1)*(2*PI)*(t-center)

def ratio(t):
 x=2*PI*t;q=expminus(x)
 U=poly(u[4:],q);V=poly(v[4:],q)
 # numerator tail bounded by ERR (same majorant as above, much larger).
 num=V-45*U/x+I(-ERR,ERR)
 # J=prod (1-q^n)^96; omitted log tail <= 192*q^41/(1-q).
 J=I(1);qp=q
 for n in range(1,41):J=J*(1-qp)**96;qp=qp*q
 rem=192*qp/(1-q)
 J=J*(1-I.raw(0,rem.h))
 return num/(990*J)

def cover(lo,hi,kind,sign,maxdepth=35):
 stack=[(Q(lo),Q(hi),0)];cells=0;bounds=[]
 while stack:
  a,b,d=stack.pop();r=station(I(a,b),kind=='D');ok=(r.l>0 if sign>0 else r.h<0)
  if ok:cells+=1;continue
  if d==maxdepth:raise RuntimeError(('unresolved',str(a),str(b),kind,sign,r.bounds()))
  m=(a+b)/2;stack.extend([(a,m,d+1),(m,b,d+1)])
 return cells

if __name__=='__main__':
 st=time.time();out={'status':'EXACT dyadic outward intervals; rational coefficient and tail bounds','precision_bits':BITS,'theta_degree':N,'tail_absolute_bound':str(ERR),'M':str(M)}
 print('tail',float(errP),float(errD),flush=True)
 # First exploratory exact point brackets.
 for ts in ['1','1.02','1.05','1.09','1.1','1.11','1.12','1.15','1.2','1.5','2']:
  print(ts,'P',station(I(ts)).bounds(),'D',station(I(ts),True).bounds(),'r',ratio(I(ts)).bounds(),flush=True)
 out['certificates']={}
 for label,a,b,kind,sign in [('near_selfdual','1','1.02','D',-1),('descending','1.02','1.09','P',-1),('minimum_window','1.09','1.12','D',1),('ascending','1.12','2.2','P',1)]:
  count=cover(a,b,kind,sign);out['certificates'][label]={'t_interval':[a,b],'quantity':kind,'sign':sign,'dyadic_subintervals':count};print(label,count,flush=True)
 lo=Q('1.09');hi=Q('1.12')
 for _ in range(55):
  md=(lo+hi)/2;p=station(I(md))
  if p.h<0:lo=md
  elif p.l>0:hi=md
  else:raise RuntimeError('root bisection uncertain')
 # For x>=x0=4.4*pi each x^j exp(-k*x), k>=1,j<=2, decreases.
 x0=2*PI*I('2.2');q0=expminus(x0)
 tailabs=x0*x0*poly([0]+[abs(c) for c in A[1:]],q0)+x0*poly([0]+[abs(c) for c in B[1:]],q0)+poly([0]+[abs(c) for c in C[1:]],q0)+I(0,ERR)
 assert A[0]==0 and B[0]==0 and C[0]*S>tailabs.h
 out['tail_monotonicity']={'starts_at_t':'2.2','positive_constant':C[0],'absolute_remainder_upper_integer':(tailabs.h+S-1)//S}
 out['unique_minimum_t_rational']=[str(lo),str(hi)];out['unique_minimum_t_display']=[float(lo),float(hi)];rr=ratio(I(lo,hi));out['minimum_r_interval']=rr.bounds()
 out['minimum_r_rational']=[str(Q(rr.l,S)),str(Q(rr.h,S))]
 scale=10**7;out['minimum_r_certified_decimal_enclosure']=[str(Q(rr.l*scale//S,scale)),str(Q(-((-rr.h*scale)//S),scale))]
 out['minimum_t_certified_decimal_enclosure']=[str(Q((lo*10**15).__floor__(),10**15)),str(Q((hi*10**15).__ceil__(),10**15))]
 out['seconds']=time.time()-st
 print(json.dumps(out,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
