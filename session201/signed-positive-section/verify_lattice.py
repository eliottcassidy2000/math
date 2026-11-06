"""Exact cubic-order and physical guard identities; Python stdlib."""
from pathlib import Path
import hashlib,json

def need(x,msg):
 if not x:raise RuntimeError(msg)
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def mul(a,b):
 c=[0]*5
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 for k in [4,3]:
  c[k-3]+=c[k];c[k-2]+=c[k];c[k-1]+=c[k]
 return tuple(c[:3])
def power(a,k):
 r=(1,0,0)
 for _ in range(k):r=mul(r,a)
 return r
def det3(A):return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
def raw(a,e,p):return a+2*e+p,e,a+2*e+2*p,e+p
def psi(s):
 D,E,U,d=s;need(U-D+E-d==0,'zero-f lattice')
 a=D-E-d;e=E;p=d-E;return a,e-p,p
def J(s):
 D,E,U,d=s;return D-E+U,D-E,D+U,D
def K(v,w):
 D,E,U,d=v;A,B,C,e=w
 return D*e+d*A+E*B+U*C-E*e-d*B-U*e-d*C+d*e
def ell(x):return x[0]+2*x[1]+4*x[2]
def main():
 theta=(0,1,0);inv=(-1,-1,1)
 need(mul(theta,inv)==(1,0,0),'integral unit inverse')
 basis=[(1,0,0),(0,1,0),(0,0,1)]
 rawbasis=[raw(1,0,0),raw(0,1,0),raw(0,1,1)]
 for x,s in zip(basis,rawbasis):
  need(psi(s)==x,'unimodular coefficient basis')
  need(psi(J(s))==mul(theta,x),'J is unit multiplication')
 gram=[[ell(mul(a,b))for b in basis]for a in basis]
 need(gram==[[1,2,4],[2,4,7],[4,7,13]]and det3(gram)==-1,'perfect pairing')
 for i in range(3):
  for j in range(3):need(K(rawbasis[i],rawbasis[j])==gram[i][j],'physical pairing equals ring pairing')
 R=[[0,0,1],[1,1,1],[0,1,0]];need(det3(R)==1,'unit determinant')
 # Characteristic polynomial checked at enough points for a monic cubic.
 for x in range(4):
  mat=[[x*int(i==j)-R[i][j]for j in range(3)]for i in range(3)]
  need(det3(mat)==x**3-x*x-x-1,'characteristic polynomial')
 need(power(theta,6)==(4,6,7)and power(theta,7)==(7,11,13),'unit powers')
 need(mul(power(theta,6),(13,-7,0))==(3,1,0),'signed boundary identity')
 z=(-1,-7,-1,-7);q=z
 for _ in range(6):q=J(q)
 need(q==(5,1,5,1)and psi(z)==(13,-7,0)and psi(q)==(3,1,0),'raw physical boundary match')
 false=(1,0,4,3);true=(1,0,4,1)
 need(J(false)==J(true)==q,'same J output different d')
 need(false[2]<2*false[3],'algebraic inverse violates physical cone')
 need(true==(1,0,4,1),'literal t.O predecessor')
 return dict(status='PASS',J_determinant=1,characteristic_polynomial='X^3-X^2-X-1',gram=gram,gram_determinant=-1,
  unit_inverse=inv,theta6=power(theta,6),theta7=power(theta,7),formal_signed_boundary=z,physical_after_six_J=q,
  nonphysical_zero_f_inverse=false,actual_off_face_predecessor=true,
  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  scope='Explicit cubic order and lattice identities; not a full ring-of-integers identification or a physical inverse without guards.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
