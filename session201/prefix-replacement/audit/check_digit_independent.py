"""Separate independent sidecar audit; original maximal-prefix audit unchanged."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib,json,subprocess,sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
FROZEN=HERE/'frozen-inputs'
EXPECTED={
 'DIGIT-ANSATZ.md':'f9a7b0537003b23b6db15423da25f8fd576461ebb58b39c210e90b037fe6fe59',
 'check_digit_ansatz.py':'c4bf47af37fd6839309d9d4def09b7ec7dee01e147c64bdd6e0979a864db56b8',
 'DIGIT-normal.json':'caeda63b8c707e6d8b4fcbe15fac125f8cca1df5e9acb2bf200938348a46ee1c',
 'DIGIT-optimized.json':'caeda63b8c707e6d8b4fcbe15fac125f8cca1df5e9acb2bf200938348a46ee1c',
 'verify.py':'2868922f90f5ba2cc8d6122408d87df692bfd64d1afa5747b814ca4f070af8dc',
}
def require(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
O=((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1))
Z=((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1))
J=((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
M={'O':O,'Z':Z,'J':J}
I=tuple(tuple(int(i==j) for j in range(4))for i in range(4))
def mm(a,b):return tuple(tuple(sum(a[i][k]*b[k][j]for k in range(4))for j in range(4))for i in range(4))
def mv(a,v):return tuple(sum(a[i][k]*v[k]for k in range(4))for i in range(4))
def matrix(w):
 out=I
 for c in w:out=mm(M[c],out)
 return out

def g(u):
 require(u>=1,'positive g index')
 answer=1
 while u>1:
  if u%2==0:u//=2
  else:
   u=(u+1)//2;answer+=u-u%2
 return answer

def op(w,u):
 for c in w:
  if c=='O':u*=2
  elif c=='Z':u=2*u-1
  elif c=='S':u+=1 if u%2 else -1
  else:raise ValueError(c)
  if u<2:return None
 return u

def words(bound):
 for l in range(bound+1):
  for w in product('OZS',repeat=l):yield ''.join(w)

def affine_from_samples(w,r):
 # Infer coefficients from independent index/g evaluations at two points.
 # This does not use the author's per-operation affine update code.
 x=100+r;xx=x+4;y=op(w,x);yy=op(w,xx)
 a=Fraction(yy-y,xx-x);b=y-a*x
 dx=g(y)-g(x);dxx=g(yy)-g(xx)
 c=Fraction(dxx-dx,xx-x);d=dx-c*x
 return a,b,c,d

def main():
 for name,h in EXPECTED.items():require(sha(FROZEN/name)==h,('frozen',name))
 seed=(1,0,0,0)
 require(mv(J,seed)==(1,1,1,1),'initial formal J')
 # Complete induction certificate: there are exactly two residues with
 # D=d=1 and U=E. Verify closure, not merely words through a depth.
 residues=[(1,0,0,1),(1,1,1,1)];closure=[]
 for s in residues:
  for name,m in M.items():
   v=tuple(x%2 for x in mv(m,s))
   require(v in residues,('complete residue closure',s,name,v))
   closure.append(dict(state=s,letter=name,output=v))
 central=0;range_ab=[None,None]
 for l in range(1,11):
  ws=['J'] if l==1 else ('J'+''.join(x)+'J'for x in product('OZJ',repeat=l-2))
  for w in ws:
   m=matrix(w);p=mm(Q,m);v=mv(m,seed)
   require(p[0][0]==m[3][0]==v[3] and p[0][0]%2==1,'AB extraction')
   require(v[0]%2==v[3]%2==1 and (v[2]-v[1])%2==0,'full integer parity')
   range_ab=[p[0][0] if range_ab[0] is None else min(range_ab[0],p[0][0]),p[0][0] if range_ab[1] is None else max(range_ab[1],p[0][0])]
   central+=1
 # Independently evaluate empty Q pairing for fixed p,q and free A,u,B,v.
 empty_checks=0
 for A,u,B,v,p,q in product(range(-2,3),repeat=6):
  s=(A,p,u,1);t=(B,q,v,1)
  actual=sum(s[i]*Q[i][j]*t[j]for i in range(4)for j in range(4))
  require(actual==A+B+u*v-u-v+p*q-p-q+1,'empty expansion')
  empty_checks+=1
 composed=author_range=0
 coeffs=[]
 for w in words(5):
  forms=[affine_from_samples(w,r)for r in range(4)]
  for u in range(2,501):
   v=op(w,u)
   if v is None:continue
   a,b,c,d=forms[u%4]
   require(v==a*u+b and g(v)==g(u)+c*u+d,('interpolated exact affine',w,u))
   composed+=1
   if len(w)<=3 and u<=300:author_range+=1
  if w in ['S','SS','SZ','ZS','OS','SO']:
   coeffs.append({'word':w,'residue_forms':[[str(x)for x in form]for form in forms]})
 require(author_range==11947,'independent author universe count')
 # General substitution invariant by matrix chain: tau' = H tau, where
 # column0 is exactly e_D and the index row has no A coefficient. Thus
 # (H_left^T P H_right)[0,0]=P[0,0] for every affine choice, exactly.
 # Verify that structural identity on all derived affine pairs for selected
 # words and all parity branches, with both original and exchanged pairing.
 matrix_controls=0
 affine_maps=[]
 for w in words(2):
  for r in range(4):
   a,b,c,d=affine_from_samples(w,r);pnew=op(w,100+r)%2
   affine_maps.append(((1,0,c,d),(0,0,0,pnew),(0,0,a,b),(0,0,0,1)))
 def tr(m):return tuple(zip(*m))
 for w in ['', 'J','JJ','JOJ','JZJ','JOZJ','JJJJ']:
  P=mm(Q,matrix(w))
  for H in affine_maps:
   for K in affine_maps:
    out=mm(mm(tr(H),P),K)
    require(out[0][0]==P[0][0],'general AB preservation control')
    require(tr(out)[0][0]==P[0][0],'exchange AB control')
    if not w:require(out[0][3]==out[3][0]==1,'empty individual A B coefficients')
    matrix_controls+=1
 for name,h in EXPECTED.items():require(sha(FROZEN/name)==h,('unchanged',name))
 return dict(status='PASS',proof_sha256=EXPECTED['DIGIT-ANSATZ.md'],verifier_sha256=sha(Path(__file__)),
  full_parity_residue_closure=closure,central_max_word_length=10,central_words=central,AB_coefficient_range=range_ab,
  empty_pairing_checks=empty_checks,affine_composition_length_bound=5,index_upper_bound=500,
  independent_affine_checks=composed,matched_author_affine_universe=author_range,
  affine_samples=coeffs,substituted_matrix_controls=matrix_controls,
  scope='Only coefficient-wise formal identities under separate index changes; no obstruction to special-index functional identities or occupancy.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
