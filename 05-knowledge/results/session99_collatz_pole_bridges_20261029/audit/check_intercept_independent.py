#!/usr/bin/env python3
"""Independent rational barriers and sibling-aware finite-height moments."""
from fractions import Fraction as F
import json

def need(p,m):
    if not p: raise RuntimeError(m)
def L(r): return F(3**r,2**(r+1)-1)-B(r)
def B(r): return F(3**r,2**(3*(r+1))-1)
barriers=[]
for r,R in [(2,F(7,5)),(3,F(19,10))]:
    f=L(r)+B(r)*R**3-R
    need(f<0,'strict root upper barrier')
    alpha=(1+R+R*R)/7
    need(alpha<1,'strict intercept contraction')
    barriers.append({'r':r,'slope_moment_upper':str(R),'alpha_upper':str(alpha),'f_at_upper':str(f)})
need(L(4)>2,'fourth moment lower barrier')

# Completed-tree truncation: q is probability all children are available.
# x=EU, m=EU², y=EV, h=EUV, t=EV², all with completion indicator.
# In particular absent siblings give factors q, not 1.
q=x=m=y=h=t=F(0)
rows=[]
for height in range(4):
    nq=F(6,7)+q**3/7
    nx=L(1)+B(1)*x**3
    nm=L(2)+B(2)*m**3
    ny=L(1)/3+x**3/63+y*(x*x+q*x+q*q)/7
    nh=L(2)/3+m**3*F(3,511)+h*(m*m+x*m+x*x)/21
    nt=L(2)/9+m**3/511+t*(m*m+q*m+q*q)/7+F(2,63)*h*(m*m+x*m+x*x)+F(2,7)*y*h*(m+x+q)
    need(all(a<=b for a,b in zip((q,x,m,y,h,t),(nq,nx,nm,ny,nh,nt))),'monotone completed trees')
    need(nq<1 and nx<1 and nm<F(7,5) and ny<F(7,12),'finite height bounds')
    # U<=3/2 on every complete tree.
    need(nm<=F(3,2)*nx and nh<=F(3,2)*ny,'pointwise bounded slope controls')
    q,x,m,y,h,t=nq,nx,nm,ny,nh,nt
    rows.append({'height':height,'completion_probability':float(q),'EU':float(x),'EU2':float(m),'EV':float(y),'EUV':float(h),'EV2':float(t)})
# Deliberately wrong omission of siblings produces a larger EV recursion.
wrong_y=L(1)/3+x**3/63+y*(x*x+x+1)/7
right_y=L(1)/3+x**3/63+y*(x*x+q*x+q*q)/7
need(wrong_y>right_y,'sibling omission must not be called equality')
print(json.dumps({'status':'PASS','barriers':barriers,'fourth_slope_moment_lower':str(L(4)),'completed_height_recursion':rows,'sibling_indicator_hostile':'strict inequality verified'},indent=2))
