"""Reconstruct sections and height Gram from primary published polynomials."""
import json,math
from pathlib import Path
import sympy as S
P=Path(__file__).resolve().parent;d=json.loads((P/'elkies_948_mathdata.json').read_text());t=S.symbols('t')
def must(x):
 if not x:raise ValueError('section audit failed')
def poly(cs):return S.Poly.from_list([S.Rational(c)for c in cs],t)
A=-27*poly(d['S']);B=S.Rational(27,4)*poly(d['T']);xs=list(map(poly,d['x_sections']));ys=[]
for x,sgn in zip(xs,d['y_leading_signs']):
 c,fs=(x**3+A*x+B).sqf_list();q=S.Rational(c);num,den=map(int,q.as_numer_denom());a,b=math.isqrt(num),math.isqrt(den);must(a*a==num and b*b==den);y=S.Poly(S.Rational(a,b),t)
 for f,e in fs:must(e%2==0);y*=f**(e//2)
 if (y.LC()>0)!=(sgn=='+'):y=-y
 must(y*y==x**3+A*x+B);ys.append(y)
must(len(xs)==len(ys)==17);must((ys[0]-poly(d['y1'])).is_zero)
D=-16*(4*A**3+27*B**2);must(D.degree()==24 and S.gcd(D,D.diff()).degree()==0)
G=S.eye(17)*4;infinity_terms=0
for i in range(17):
 for j in range(i):
  dx=xs[i]-xs[j];dy=ys[i]-ys[j];inf=min(4-dx.degree(),6-dy.degree());infinity_terms+=int(inf)
  h=2-S.gcd(dx,dy).degree()-inf;G[i,j]=G[j,i]=h
must(G==S.Matrix(d['published_gram']));must(G.det()==948)
minors=[int(G[:i,:i].det())for i in range(1,18)];must(all(x>0 for x in minors))
for k,m in enumerate(d['secant_slopes'],1):must((ys[k]-ys[0]-poly(m)*(xs[k]-xs[0])).is_zero)
print(json.dumps({'sections':17,'gram_matches':True,'determinant':int(G.det()),'leading_principal_minors':minors,'homogeneous_infinity_intersections':infinity_terms,'verified_quadratic_slopes':15},indent=2))
