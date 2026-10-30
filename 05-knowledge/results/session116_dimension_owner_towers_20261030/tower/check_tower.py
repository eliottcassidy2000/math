"""Exact CD doubling, polynomial defect and scoped zero-product checks; stdlib."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
COUNT=0
def require(ok,msg):
 global COUNT
 COUNT+=1
 if not ok: raise RuntimeError(msg)
def add(x,y): return tuple(a+b for a,b in zip(x,y))
def neg(x): return tuple(-a for a in x)
def sub(x,y): return add(x,neg(y))
def scale(k,x):return tuple(k*a for a in x)
def conj(x):return (x[0],)+tuple(-a for a in x[1:])
def norm(x):return sum(a*a for a in x)
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def mul(x,y):
 n=len(x)
 if n==1:return (x[0]*y[0],)
 h=n//2;a,b=x[:h],x[h:];c,d=y[:h],y[h:]
 return sub(mul(a,c),mul(conj(d),b))+add(mul(d,a),mul(b,conj(c)))
def basis(n,i):return tuple(int(j==i) for j in range(n))
def assoc(a,b,c):return sub(mul(mul(a,b),c),mul(a,mul(b,c)))
def W(a,s):return a+(0,)*8+s

def main():
 # Quartic cross term coefficients: independently evaluate the two mixed inner products.
 for i,j,k,l in product(range(8),repeat=4):
  a,b,c,d=(basis(8,z) for z in (i,j,k,l))
  cross=2*(dot(mul(d,a),mul(b,conj(c)))-dot(mul(a,c),mul(conj(d),b)))
  require(cross==2*dot(assoc(d,a,c),b),'associator coefficient')
 # A direct exact norm polynomial check for all 2-term signed input vectors,
 # including repeated units. It also verifies the dimension boundary 8 vs16.
 for n in (4,8,16):
  vec=[add(basis(n,i),scale(s,basis(n,j))) for i in range(n) for j in range(i,n) for s in (-1,1)]
  defects=set()
  for x in vec:
   for y in vec:
    delta=norm(mul(x,y))-norm(x)*norm(y);defects.add(delta)
    if n<=8:require(delta==0,'composition below16')
    else:
     a,b=x[:8],x[8:];c,d=y[:8],y[8:]
     require(delta==2*dot(assoc(d,a,c),b),'full defect')
  require(defects==({0} if n<=8 else {-4,0,4}),'sparse defect boundary')
 # All quaternion factors from the 9-element signed basis-plus-zero set.
 bank=[(0,)*4]+[scale(s,basis(4,i)) for i in range(4) for s in (-1,1)]
 zeros=0
 for a,s,c,t in product(bank,repeat=4):
  x,y=W(a,s),W(c,t)
  actual=mul(x,y)==(0,)*16
  predicted=(norm(x)==0 or norm(y)==0)
  if norm(x) and norm(y) and norm(a):
   predicted=(a[0]==c[0]==0 and dot(a,c)==0 and norm(a)==norm(s) and norm(c)>0 and scale(norm(a),t)==mul(mul(s,a),c))
  require(actual==predicted,'zero-product classification')
  expected=sub(mul(a,c),mul(conj(s),t))+(0,)*8+add(mul(t,conj(a)),mul(s,c))
  require(mul(x,y)==expected,'quaternion reduction')
  if actual and norm(x) and norm(y):zeros+=1
 # Integer cone samples, including arbitrary norm rather than only unit a.
 integer_families=0
 for A in product(range(-2,3),repeat=3):
  if A==(0,0,0):continue
  a=(0,)+A;s=a;N=norm(a)
  for C in product(range(-2,3),repeat=3):
   if C==(0,0,0) or sum(u*v for u,v in zip(A,C)):continue
   c=(0,)+C;x=W(a,s);y=W(scale(N,c),mul(mul(s,a),c))
   require(mul(x,y)==(0,)*16,'integer cone family')
   inverse=scale(F(1,norm(x)),conj(x))
   require(mul(inverse,x)==basis(16,0) and mul(x,inverse)==basis(16,0),'two-sided inverse')
   require(assoc(inverse,x,y)==y,'cancellation associator')
   integer_families+=1
 # Same explicit pair embeds unchanged in every subsequent CD stage tested.
 i,j,k=(basis(4,z) for z in (1,2,3));one=basis(4,0)
 x,y=W(i,one),W(j,k)
 require(mul(x,y)==(0,)*16,'named pair')
 for n in (16,32,64,128):
  xx=x+(0,)*(n-16);yy=y+(0,)*(n-16)
  require(mul(xx,yy)==(0,)*n and norm(xx)==norm(yy)==2,'embedding')
  require(mul(xx,conj(xx))==scale(norm(xx),basis(n,0)),'conjugation norm')
 out={'status':'PASS','gates':COUNT,'signed_basis_nontrivial_zero_products':zeros,'integer_cone_pairs':integer_families,'embedding_dimensions':[16,32,64,128],'sparse_defects':[-4,0,4],'proof_scope':'all-real theorem proved in TOWER_DEFECT.md; tests finite-exact'}
 Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
