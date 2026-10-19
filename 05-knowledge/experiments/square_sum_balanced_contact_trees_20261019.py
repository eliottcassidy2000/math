#!/usr/bin/env python3
"""Exact fixed-map, free-jet, analytic-radius, and graph controls."""
from fractions import Fraction as F
from math import isqrt,gcd,lcm,ceil
from collections import Counter
from hashlib import sha256
import json,sys
import sympy as s
sys.set_int_max_str_digits(0)
checks=Counter()
def check(v,n):
 checks[n]+=1
 if not v:raise RuntimeError(n)
def add(R,S,k):
 if R is None:return S
 if S is None:return R
 x,y=R;u,v=S
 if x==u and y==-v:return None
 slope=(3*x*x-2*k*k*x-1)/(2*y) if R==S else (v-y)/(u-x)
 xx=slope*slope+k*k-x-u
 return xx,-y+slope*(x-xx)
def mul(R,n,k):
 out=None
 while n:
  if n&1:out=add(out,R,k)
  R=add(R,R,k);n//=2
 return out
def anchor(t):return ((t*t+4*t+1)/(t*t+1),8*t*(t+1)/(t*t+1)**2)
def kval(t):return (t*t+2*t-1)/(t*t+1)
def f(t):return 4*t*(t*t-1)/(t*t+1)**2
def parameter(R,k):
 x,y=R
 return (k*y+x-k*k)/(x*(x-k*k))
t,h=s.symbols('t h');k=kval(t);Q=anchor(t);P=(s.S(0),k)
def sadd(V,W):return tuple(s.factor(s.cancel(z)) for z in add(V,W,k))
V=sadd(Q,P);S=sadd(V,V);R3=sadd(Q,S);R5=sadd(R3,S)
F3=s.factor(s.cancel(parameter(R3,k)));F5=s.factor(s.cancel(parameter(R5,k)))
functions={(1,3):F3,(1,5):F5,(-1,3):s.factor(s.cancel(1/F3.subs(t,1/t))),(-1,5):s.factor(s.cancel(1/F5.subs(t,1/t)))}
expected={(1,3):[1,3,0,2],(1,5):[1,5,0,0],(-1,3):[1,3,6,14],(-1,5):[1,5,20,80]}
mapdata={};polys=[];r0=F(1,64);K=6
for key,fn in functions.items():
 N,D=s.fraction(s.cancel(fn));degree=max(s.degree(N,t),s.degree(D,t))
 check(degree==(11 if key[1]==3 else 33),'fixed_map_degrees')
 check(s.degree(s.gcd(N,D),t)==0,'reduced_map_forms')
 nn,dd=s.Poly(N.subs(t,1+h),h),s.Poly(D.subs(t,1+h),h)
 if dd.nth(0)<0:nn,dd=-nn,-dd
 co=[]
 for j in range(4):co.append(s.cancel((nn.nth(j)-sum(dd.nth(i)*co[j-i] for i in range(1,j+1)))/dd.nth(0)))
 check(co==expected[key],'exact_four_term_jet')
 lower,rem=s.div(nn-(1+h)*dd,s.Poly(h,h));check(rem.is_zero,'cusp_linear_factor')
 U,rem=s.div(nn-dd,s.Poly(h,h));check(rem.is_zero,'displacement_linear_factor')
 tests=[dd,lower,5*dd-4*nn,nn.diff()*dd-nn*dd.diff()]
 for p in tests:
  p0=F(p.nth(0));tail=sum(F(abs(p.nth(i))) for i in range(1,p.degree()+1))
  check(p0>0,'positive_chart_constant')
  if tail:r0=min(r0,p0/(2*tail))
  polys.append(p)
 K=max(K,ceil(2*sum(F(abs(z)) for z in U.all_coeffs())/F(dd.nth(0))))
 mapdata[str(key)]={'degree':int(degree),'jet':[str(c) for c in co],
     'numerator_coefficients':[int(z) for z in s.Poly(N,t).all_coeffs()],
     'denominator_coefficients':[int(z) for z in s.Poly(D,t).all_coeffs()]}
for p in polys:
 tail=sum(F(abs(p.nth(i)))*r0**i for i in range(1,p.degree()+1))
 check(tail<=F(p.nth(0))/2,'uniform_real_and_complex_cutoff')
# Prime-valuation comparison for the d=5 section.
N5,D5=s.fraction(s.cancel(F5));rest,rem=s.div(N5*N5+D5*D5,t*t+1,t)
check(rem==0,'F5_hypotenuse_factor')
check(s.rem(rest,t*t+1,t)==-2**52,'F5_exact_old_prime_preservation')
check(s.resultant(N5,D5,t)==2**968,'F5_odd_prime_good_reduction')
def compose(a,b):
 A,B,C=a;x,y,z=b
 return A*x,A*y+B*x*x,A*z+2*B*x*y+C*x*x*x
jet_records=[]
for eps in [1,-1]:
 f3=tuple(expected[(eps,3)][1:]);f5=tuple(expected[(eps,5)][1:])
 A=compose(f3,f5);B=compose(f5,f3)
 check(A==((15,0,250) if eps==1 else (15,210,3190)),'block_A_jet')
 check(B==((15,0,10) if eps==1 else (15,210,2950)),'block_B_jet')
 nodes={'':(1,0,0)};level={'':(1,0,0)}
 for depth in range(6):
  new={}
  for word,j in level.items():
   for letter,first,block in [('A',f5,A),('B',f3,B)]:
    nodes[word+letter+'*']=compose(first,j)
    new[word+letter]=compose(block,j)
    if eps==1:
     digits=word+letter
     expected_c=15**(len(digits)-1)*sum((250 if c=='A' else 10)*225**i for i,c in enumerate(digits))
     check(new[word+letter][2]==expected_c,'base225_formula_vs_composition')
  nodes.update(new);level=new
 check(len(set(nodes.values()))==len(nodes),'all_prefix_jets_distinct')
 vals=list(nodes.values())
 for i,u in enumerate(vals):
  for v in vals[i+1:]:
   check(u[0]!=v[0] or (u[1]==v[1] and u[2]!=v[2]),'first_distinction_order_one_or_three')
 jet_records.append({'epsilon':eps,'block_depth':6,'prefix_count':len(nodes),'all_jets_distinct':True})
def eval_poly(coeffs,x):
 out=F(0)
 for z in coeffs:out=out*x+z
 return out
def eval_map(eps,d,x):
 rec=mapdata[str((eps,d))]
 return eval_poly(rec['numerator_coefficients'],x)/eval_poly(rec['denominator_coefficients'],x)
def square(r):
 a,b=isqrt(r.numerator),isqrt(r.denominator)
 check(a*a==r.numerator and b*b==r.denominator,'rational_contact_square')
 return F(a,b)
for eps in [1,-1]:
 for d in [3,5]:
  for x in [1+F(1,100),1+F(1,1000)]:
   y=eval_map(eps,d,x)
   check(x<y and 0<f(x)<f(y)<F(1,2),'positive_finite_section_control')
   square(1+eps*f(x)-eps*f(y))
   # Independent rational group expression dQ+(d-1)P on the plus input,
   # including reciprocal input for the minus formula.
   xx=x if eps==1 else 1/x;kk=kval(xx)
   R=add(mul(anchor(xx),d,kk),mul((F(0),kk),d-1,kk),kk)
   u=parameter(R,kk);u=u if eps==1 else 1/u
   check(u==y,'independent_rational_section_group_law')
# Explicit once-subdivided depth-one binary contact tree, both ports.
graphs=[]
for eps in [1,-1]:
 root=1+F(1,10**8);v3=eval_map(eps,3,root);v5=eval_map(eps,5,root)
 A=eval_map(eps,3,v5);B=eval_map(eps,5,v3)
 ts=[root,v5,A,v3,B];edges=[(0,1),(1,2),(0,3),(3,4)]
 check(len(set(ts))==5,'numeric_tree_parameters_distinct')
 cs=[];triples=[]
 for x in ts:
  n,m=x.numerator,x.denominator;delta=2 if n%2 and m%2 else 1
  a,b,c=(n*n-m*m)//delta,2*n*m//delta,(n*n+m*m)//delta
  check(a*a+b*b==c*c and gcd(a,b)==1 and c*c>4*a*b,'numeric_tree_primitive_admissible')
  triples.append((a,b,c));cs.append(c)
 H=lcm(*cs);M=2*H*H;ports=[]
 for a,b,c in triples:
  g=H//c;ports.append((2*g*g*(c*c-4*a*b),2*g*g*(c*c+4*a*b)))
 labels=[M]+[x for pair in ports for x in pair]
 check(len(set(labels))==11 and min(labels)>0,'integer_tree_labels_distinct_positive')
 integer_edges=[]
 for lo,hi in ports:integer_edges.extend([(lo,M),(lo,hi),(M,hi)])
 for i,j in edges:
  x=ports[i][1 if eps==1 else 0];y=ports[j][0 if eps==1 else 1]
  integer_edges.append((x,y));square(1+eps*f(ts[i])-eps*f(ts[j]))
 for a,b in integer_edges:check(isqrt(a+b)**2==a+b,'integer_tree_square_edge')
 encoded=json.dumps({'labels':labels,'edges':integer_edges},separators=(',',':'))
 graphs.append({'epsilon':eps,'petals':5,'integer_vertices':11,'checked_edges':len(integer_edges),
   'hypotenuse_bits':[c.bit_length() for c in cs],'largest_label_bits':max(labels).bit_length(),
   'certificate_sha256':sha256(encoded.encode()).hexdigest()})
print(json.dumps({'status':'FINITE-EXACT; all-size theorem proved separately','checks':dict(checks),
 'analytic_constants':{'r0':str(r0),'K':K},'maps':mapdata,'jet_tree_controls':jet_records,
 'integer_tree_controls':graphs},indent=2))
