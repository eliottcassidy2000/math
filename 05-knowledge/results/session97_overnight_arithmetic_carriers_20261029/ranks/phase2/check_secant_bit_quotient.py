"""15 secant ambiguity bits: free for basis rank, observable for named sections."""
import json,math
from pathlib import Path
import sympy as S
P=Path(__file__).resolve().parent;d=json.loads((P/'elkies_948_mathdata.json').read_text());t=S.symbols('t')
def must(v):
 if not v:raise ValueError('secant quotient audit failed')
def poly(cs):return S.Poly.from_list([S.Rational(x)for x in cs],t)
A=-27*poly(d['S']);B=S.Rational(27,4)*poly(d['T']);xs=list(map(poly,d['x_sections']));x1=xs[0];y1=poly(d['y1']);ys=[y1];alts=[]
for x,m in zip(xs[1:16],map(poly,d['secant_slopes'])):
 y=y1+m*(x-x1);rx=m*m-x1-x;ry=y1+m*(rx-x1)
 must((y*y-x**3-A*x-B).is_zero);must((ry*ry-rx**3-A*rx-B).is_zero)
 must(not (rx-x).is_zero);must((ry-y1-m*(rx-x1)).is_zero);alts.append((rx,ry));ys.append(y)
# Remaining section is reconstructed solely for named-pairing controls.
c,fs=(xs[16]**3+A*xs[16]+B).sqf_list();q=S.Rational(c);num,den=map(int,q.as_numer_denom());y=S.Poly(S.Rational(math.isqrt(num),math.isqrt(den)),t)
for f,e in fs:must(e%2==0);y*=f**(e//2)
if (y.LC()>0)!=(d['y_leading_signs'][16]=='+'):y=-y
must((y*y-xs[16]**3-A*xs[16]-B).is_zero);ys.append(y)
G=S.Matrix(d['published_gram']);masks=[0,(1<<15)-1,0x5555,0x2aaa]+[1<<k for k in range(15)];changed=0
for mask in masks:
 U=S.eye(17);xx=list(xs);yy=list(ys)
 for k in range(15):
  if mask>>k&1:U[k+1,k+1]=-1;U[0,k+1]=-1;xx[k+1],yy[k+1]=alts[k]
 H=U.T*G*U;must(U.det()==(-1)**mask.bit_count());must(H.det()==948)
 for i in range(17):
  for j in range(i):
   dx=xx[i]-xx[j];dy=yy[i]-yy[j];h=2-S.gcd(dx,dy).degree()-min(4-dx.degree(),6-dy.degree());must(h==H[i,j])
 if H!=G:changed+=1
print(json.dumps({'ambiguous_generators':15,'distinct_named_bases_by_unimodular_argument':2**15,'masks_with_full_height_replay':len(masks),'changed_gram_matrices':changed,'all_determinants':948},indent=2))
