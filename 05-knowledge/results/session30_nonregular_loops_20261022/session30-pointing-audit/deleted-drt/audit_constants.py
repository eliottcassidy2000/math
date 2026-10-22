"""Independent leading contraction and direct differential-operator check."""
import json, math
from pathlib import Path
from collections import Counter
import sympy as s
base=Path(__file__).resolve().parent
source=base.parents[1]/'session30-decorated-engine/symbolic-kernels.json'
z,u=s.symbols('z u')
cores=json.loads(source.read_text());answers=[];A2=0
for c in cores:
 if c['weight']!=2:continue
 es=c['edges'];degree=[0]*c['V']
 for a,b in es:degree[a]+=1;degree[b]+=1
 terminal=[d==1 for d in degree]
 assert len(es)==c['V']+2-sum(terminal)
 internal=Counter();colors=[0]*c['V'];factors=[]
 for a,b in es:
  assert not(terminal[a] and terminal[b])
  if terminal[a] or terminal[b]:
   sign=-1 if terminal[a] else 1
   factors.append([sign*z,-sign*z,0,0]);colors[b if terminal[a] else a]+=1
  elif a==b:factors.append([-z*z,-z*z,z*z,-z*z])
  else:
   assert a<b
   internal[(a,b)]+=1;factors.append([z,-z,0,0])
 survived=not(any(m%2 for m in internal.values()) or any(m%2 for m in colors))
 val=0
 if survived:
  for term in c['terms']:
   val+=s.Rational(term['numerator'],c['denominator'])*s.prod(factors[i][j] for i,j in enumerate(term['indices']))
 val=s.expand(val);A2+=val
 if val:answers.append({'edges':es,'value':str(val)})
assert s.expand(A2)==7*z**4+26*z**6
# Evaluate Euler polynomial actions directly, without constructing Bell polynomials.
D=lambda f:z*s.diff(f,z)
X0=z*z*(1-u)+z**6*u*u/3
exp0=s.exp(X0)
der=[exp0]
for _ in range(4):der.append(D(der[-1]))
p1=s.simplify((der[2]-der[1])/exp0/2)
p2=s.simplify((3*der[4]-2*der[3]-3*der[2]+2*der[1])/exp0/24)
A1=-(z*z+s.Rational(7,2)*z**4)*u+(12*z**4-s.Rational(8,3)*z**6)*u*u
# Product Euler action builds the mixed term, not the producer's 'shift' formula.
mixed=s.simplify((D(D(exp0*A1))-D(exp0*A1))/exp0/2)
truncate=lambda x:s.series(x.subs(z,1),u,0,3).removeO().expand()
Y1=truncate(A1+u*p1)
Y2=truncate(u*u*A2+A1*A1/2+u*mixed+u*u*p2)
K1=truncate(s.exp(X0-1)*(1+Y1))
K2=truncate(s.exp(X0-1)*(1+Y1+Y2))
assert Y1==-s.Rational(3,2)*u+s.Rational(13,3)*u*u
assert Y2==s.Rational(35,24)*u*u
assert K1==1-s.Rational(5,2)*u+s.Rational(20,3)*u*u
assert K2==1-s.Rational(5,2)*u+s.Rational(65,8)*u*u
out={'status':'PASS','survivors':answers,'A2_coefficient':str(A2),'Y1':str(Y1),'Y2':str(Y2),'kappa1_over_e':str(K1),'kappa2_over_e':str(K2),'scope':'finite five/26-core predictors only; exact Hamiltonian asymptotics not claimed'}
(base/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
