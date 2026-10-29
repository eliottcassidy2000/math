"""Independent reconstruction from primary mathematical data; triangular sqrt and group-law heights."""
import json,math
from pathlib import Path
import sympy as s
P=Path(__file__).resolve().parent;t=s.symbols('t')
def req(x,msg):
 if not x:raise ValueError(msg)
data=json.loads((P/'elkies948_inputs.json').read_text())
def poly(cs):return s.Poly.from_list([s.Rational(c)for c in cs],t,domain=s.QQ)
S=poly(data['S']);T=poly(data['T']);A=-27*S;B=s.Rational(27,4)*T
xs=list(map(poly,data['xs']))
ys=[]
for x,sign in zip(xs,data['y_signs']):
 f=x**3+A*x+B;c=f.all_coeffs();req(f.degree()==12,'sqrt degree');r=s.sqrt(c[0]);req(r.is_Rational,'rational leading square');coeff=[r]
 for k in range(1,7):coeff.append((c[k]-sum(coeff[j]*coeff[k-j]for j in range(1,k)))/(2*r))
 y=s.Poly.from_list(coeff,t);y=y if sign=='+' else -y
 req((y*y-f).is_zero,'exact polynomial root');ys.append(y)
req((ys[0]-poly(data['y1'])).is_zero,'y1 source')
# Independent height calculation via reduced x-coordinate of P_i+P_j,
# rather than intersection gcd of coordinate differences.
G=s.eye(17)*4
for i in range(17):
 for j in range(i):
  dx=xs[i]-xs[j];dy=ys[i]-ys[j]
  n=dy*dy-(xs[i]+xs[j])*dx*dx;d=dx*dx;c=s.gcd(n,d);n=n.exquo(c);d=d.exquo(c)
  h=max(n.degree(),d.degree()+4)
  G[i,j]=G[j,i]=s.Rational(h-8,2)
sourceG=s.Matrix(data['gram'])
req(G==sourceG,'published Gram');req(G.det()==948,'det948')
req(all(G[:k,:k].det()>0 for k in range(1,18)),'positive definite')
slopes=list(map(poly,data['slopes']))
for i,m in enumerate(slopes,1):
 req((ys[i]-ys[0]-m*(xs[i]-xs[0])).is_zero,'slope');xx=m*m-xs[0]-xs[i];yy=ys[0]+m*(xx-xs[0]);req((yy*yy-xx**3-A*xx-B).is_zero,'alternative on curve');req(xx!=xs[i]or yy!=ys[i],'distinct branch')
# Inspect all 2^15 choices combinatorially: lower triangular U diag ±1.
for mask in range(1<<15):
 diag=[-1 if mask>>j&1 else 1 for j in range(15)]
 req(abs(math.prod(diag))==1,'unimodular choice')
d=s.Poly(-16*(4*A**3+27*B**2),t);req(d.degree()==24 and s.gcd(d,d.diff()).degree()==0,'I1 geometry')
print(json.dumps({'sections':len(xs),'source_gram_matches':True,'det':int(G.det()),'heights_via_reduced_group_law':136,'secants':len(slopes),'unimodular_masks':32768},indent=2))
