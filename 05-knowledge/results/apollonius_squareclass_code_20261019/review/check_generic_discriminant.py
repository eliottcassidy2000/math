#!/usr/bin/env python3
"""Independent direct-coordinate Apollonius review; no coordinate rotation."""
import sympy as s, itertools,json
from fractions import Fraction

def require(condition,message):
 if not condition:raise RuntimeError(message)
a,b,c,d,r,u,v=s.symbols('a b c d r u v')
area=a*d-b*c
q1=(a*a+b*b+r*r-u*u)/2;q2=(c*c+d*d+r*r-v*v)/2
t1=r-u;t2=r-v
z0=s.Matrix([d*q1-b*q2,-c*q1+a*q2]);z1=s.Matrix([d*t1-b*t2,-c*t1+a*t2])
# Numerator of f(R) after multiplication by area^2.
A=(z1.dot(z1)-area**2);B=2*(z0.dot(z1)-r*area**2);C=z0.dot(z0)-r*r*area**2
prod=(a*a+b*b-(r-u)**2)*(c*c+d*d-(r-v)**2)*((a-c)**2+(b-d)**2-(u-v)**2)
require(s.Poly(s.expand(B*B-4*A*C- prod*area**2),a,b,c,d,r,u,v).is_zero, "general-coordinate discriminant identity")
# Therefore disc f = product / area^2.
signs=[(1,x,y) for x,y in itertools.product((1,-1),repeat=2)]
pairs=list(itertools.combinations(range(3),2));rows=[]
for e in signs:
 row=sum(1<<(2*j+(e[i]*e[k]==-1)) for j,(i,k) in enumerate(pairs));rows.append(row)
rels=[mask for mask in range(16) if __import__('functools').reduce(int.__xor__,(rows[i] for i in range(4) if mask>>i&1),0)==0]
require(rels==[0,15], "unique squareclass relation")
# Obtain image from all six independent defect square roots, without rank oracle.
images={tuple(sum(((rows[i]>>j)&1)*((eps>>j)&1) for j in range(6))%2 for i in range(4)) for eps in range(64)}
require(images=={e for e in itertools.product((0,1),repeat=4) if sum(e)%2==0}, "even sign code")
counts={k:sum(2*(4-sum(e))==k for e in images) for k in (0,4,8)}
require(counts=={0:1,4:6,8:1}, "fixed-point histogram")
print(json.dumps({'status':'PASS','general_center_discriminant':'product defects / determinant(center differences)^2','normalization':'No rotation or square-root base extension used','sign_rows':rows,'relations':rels,'group_order':len(images),'fixed_point_counts':counts,'review_boundary':'Generic over labelled circles with radii in the base field; no numerical Frobenius distribution conclusion'},indent=2))
