#!/usr/bin/env python3
"""Exact periodic-tent and d-block primitives for the cell-90 referee."""
from bisect import bisect_right
from fractions import Fraction as F

L=168
CELL=90

def require(condition,detail):
 if not condition:raise RuntimeError(detail)

def tent_value(x,z,w):
 """Normalized pair-overlap tent at circular displacement x."""
 x=F(x)%1;u=min(x,1-x)
 outer=F(w+z,14*z);inner=F(w-z,14*z)
 height=F(24,w);slope=F(168,w)
 if u<=inner:return height
 if u>=outer:return F(0)
 return slope*(outer-u)

def tent_slope(x,z,w):
 x=F(x)%1;outer=F(w+z,14*z);inner=F(w-z,14*z)
 slope=F(168,w)
 if inner<x<outer:return -slope
 if 1-outer<x<1-inner:return slope
 return F(0)

class PeriodicPL:
 """Exact piecewise-linear G_d(x)=sum_{r<d} f(x+r alpha)."""
 def __init__(self,z,w,alpha,d):
  self.z=F(z);self.w=F(w);self.alpha=F(alpha)%1;self.d=d
  outer=F(w+z,14*z);inner=F(w-z,14*z)
  raw={F(0)}
  for r in range(d):
   shift=r*self.alpha
   for c in (inner,outer,1-outer,1-inner):raw.add((c-shift)%1)
  self.bp=sorted(raw)+[F(1)]
  self.value=[];self.slope=[];self.prefix_values=[F(0)]
  for x,y in zip(self.bp,self.bp[1:]):
   mid=(x+y)/2
   value=sum((tent_value(x+r*self.alpha,z,w) for r in range(d)),F(0))
   slope=sum((tent_slope(mid+r*self.alpha,z,w) for r in range(d)),F(0))
   self.value.append(value);self.slope.append(slope)
   width=y-x
   self.prefix_values.append(self.prefix_values[-1]+value*width+slope*width*width/2)
  self.mean=self.prefix_values[-1]
  require(self.mean==d*F(24,7*z),("mean",z,w,alpha,d,self.mean))
  values=tuple(self.at(x) for x in self.bp[:-1])
  self.oscillation=max(values)-min(values)
  self.derivative_variation=sum(abs(self.slope[i]-self.slope[i-1])
                                for i in range(len(self.slope)))

 def at(self,x):
  x=F(x)%1;i=bisect_right(self.bp,x)-1;i=min(i,len(self.value)-1)
  return self.value[i]+self.slope[i]*(x-self.bp[i])

 def prefix(self,x):
  x=F(x);require(0<=x<=1,("prefix domain",x))
  if x==1:return self.prefix_values[-1]
  i=bisect_right(self.bp,x)-1;i=min(i,len(self.value)-1)
  width=x-self.bp[i]
  return self.prefix_values[i]+self.value[i]*width+self.slope[i]*width*width/2

 def interval(self,start,length):
  start=F(start)%1;length=F(length)
  require(0<=length<=1,("interval length",length))
  if start+length<=1:return self.prefix(start+length)-self.prefix(start)
  return self.mean-self.prefix(start)+self.prefix(start+length-1)

def block_data(e,p,f,q,d):
 require(p<q<2*p,("cap-two levels",e,p,f,q,d))
 z=L*p-e;w=L*q-f;alpha=F(w-z,z);n=p//d
 beta=(d*alpha)%1
 require(beta!=0,("zero d-step",e,p,f,q,d))
 forward=beta<=F(1,2);rho=beta if forward else 1-beta
 R=e*CELL%L;S=f*CELL%L;det=R*w-S*z
 require(det%L==0,("nonintegral phase",e,p,f,q,d))
 x0=F(det//L,z)%1
 return z,w,alpha,n,forward,rho,x0

def actual_phase_lower(e,p,f,q,d):
 """Rigorous composite-trapezoid lower bound for d*floor(p/d) samples."""
 z,w,alpha,n,forward,rho,x0=block_data(e,p,f,q,d)
 G=PeriodicPL(z,w,alpha,d);T=n*rho;m=T.numerator//T.denominator;theta=T-m
 start=x0 if forward else x0-theta
 integral=m*G.mean+G.interval(start,theta)
 xn=x0+theta if forward else x0-theta
 endpoint=(G.at(x0)-G.at(xn))/2
 # Each derivative-kink phase has at most ceil(T) lifts on the open path.
 ceilT=-((-T.numerator)//T.denominator)
 error=rho*ceilT*G.derivative_variation/8
 return integral/rho+endpoint-error

def direct_complete_sum(e,p,f,q,d):
 """Literal sum of the same complete d-block subset; hostile controls only."""
 z,w,alpha,n,forward,rho,x0=block_data(e,p,f,q,d)
 G=PeriodicPL(z,w,alpha,d);sign=1 if forward else -1
 return sum((G.at(x0+sign*k*rho) for k in range(n)),F(0))

def peano_complete_sum(e,p,f,q,d):
 """Exact signed Peano identity, independent of direct sample iteration."""
 z,w,alpha,n,forward,rho,x0=block_data(e,p,f,q,d)
 G=PeriodicPL(z,w,alpha,d);T=n*rho;m=T.numerator//T.denominator;theta=T-m
 start=x0 if forward else x0-theta
 result=(m*G.mean+G.interval(start,theta))/rho
 xn=x0+theta if forward else x0-theta
 result+=(G.at(x0)-G.at(xn))/2
 correction=F(0)
 for i,bp in enumerate(G.bp[:-1]):
  jump=G.slope[i]-G.slope[i-1]
  if jump==0:continue
  first=((bp-x0) if forward else (x0-bp))%1
  lift=0 if first else 1
  while first+lift<T:
   t=first+lift;k=t//rho;u=(t-k*rho)/rho
   correction+=jump*rho*u*(1-u)/2
   lift+=1
 return result+correction
