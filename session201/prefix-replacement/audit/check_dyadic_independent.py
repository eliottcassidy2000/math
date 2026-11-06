"""Independent exact audit of actual dyadic-class rank and product rank."""
from itertools import product
from pathlib import Path
import hashlib,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
FROZEN=HERE/'frozen-inputs'
EXPECTED={
 'DYADIC-SPAN.md':'efb7ad8eefe15b68bd3eb7862ec585b03ef88f71a11a36822c2d83bc5bf4fbc7',
 'check_dyadic_span.py':'10f9902931ad1568df127c3313c019fb028bd6484dba7719d0f04591bb4fd9e0',
 'DYADIC-normal.json':'413b7fcc8f0890dea35a911ded655ed584e5496f8da975fbf7112d9ff2ca7b32',
 'DYADIC-optimized.json':'413b7fcc8f0890dea35a911ded655ed584e5496f8da975fbf7112d9ff2ca7b32',
}
def need(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def g(u):
 value=1
 while u>1:
  if u%2==0:u//=2
  else:u=(u+1)//2;value+=u-u%2
 return value

def coefficients(bits):
 # Closed sums, independent of the author's sequential coefficient update.
 k=len(bits)
 A=2**k
 r=sum(bit*2**(k-1-i) for i,bit in enumerate(bits))
 a=sum(bit*2**i for i,bit in enumerate(bits))
 b=-bits[0] if bits else 0
 c=-sum(bits[i]*bits[j]*2**(i-1-j) for i in range(k) for j in range(i))
 c-=sum(bits[i]*bits[i-1] for i in range(1,k))
 return A,r,a,b,c

def image(bits,h):
 u=h
 for bit in bits:u=2*u-bit
 return u

def det3(a):
 return sum((-1 if (i,j,k) in [(0,2,1),(1,0,2),(2,1,0)] else 1)*a[0][i]*a[1][j]*a[2][k]
            for i,j,k in [(0,1,2),(1,2,0),(2,0,1),(0,2,1),(1,0,2),(2,1,0)])

def det_bareiss(a):
 a=[row[:] for row in a];n=len(a);den=1;sign=1
 for k in range(n-1):
  if a[k][k]==0:
   row=next((i for i in range(k+1,n)if a[i][k]),None)
   if row is None:return 0
   a[k],a[row]=a[row],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    value=a[i][j]*pivot-a[i][k]*a[k][j]
    need(value%den==0,'fraction-free determinant division')
    a[i][j]=value//den
   a[i][k]=0
  den=pivot
 return sign*a[-1][-1]

def main():
 for name,h in EXPECTED.items():need(sha(FROZEN/name)==h,('frozen',name))
 blocks=profiles=spans=0;small_bases=[]
 for k in range(11):
  residues=set()
  for bits in product([0,1],repeat=k):
   A,r,a,b,c=coefficients(bits)
   need(b in (0,-1)and 0<=r<A,'coefficient structure');residues.add((-r)%A);blocks+=1
   for h in range(2,101):
    u=image(bits,h)
    need(u==A*h-r and g(u)==g(h)+a*h+b*(h%2)+c,'closed sum block formula');profiles+=1
   for m in range(1,21):
    hs=[2**m,2**(m+1),2**(m+2)-1]
    need([g(h)for h in hs]==[1,1,2**(m+1)+1],'special actual g values')
    base=[(g(h)+b*(h%2),h,1)for h in hs]
    points=[(g(image(bits,h)),image(bits,h),1)for h in hs]
    exact=-2**(k+m)*(2**(m+1)+b)
    need(det3(points)==exact==A*det3(base)and exact!=0,'all-order determinant formula');spans+=1
    if k<=3 and m==3:small_bases.append({'k':k,'residue':(-r)%A,'points':points,'det':exact})
  need(residues==set(range(A)),'every residue represented exactly')
 # Full 9-by-9 evaluation matrices for independently varying index pairs.
 # This checks the exact product inference, not just separate spans.
 products=0
 for left in small_bases:
  for right in small_bases:
   evaluation=[[x[i]*y[j]for i in range(3)for j in range(3)]
               for x in left['points']for y in right['points']]
   expected=left['det']**3*right['det']**3
   need(det_bareiss(evaluation)==expected!=0,'independent product full rank')
   products+=1
 for name,h in EXPECTED.items():need(sha(FROZEN/name)==h,('unchanged',name))
 return dict(status='PASS',proof_sha256=EXPECTED['DYADIC-SPAN.md'],verifier_sha256=sha(Path(__file__)),
  binary_block_length_bound=10,binary_blocks=blocks,profile_seed_range=[2,100],profile_controls=profiles,
  power_parameter_range=[1,20],exact_span_determinants=spans,
  determinant_formula='-2^(k+m)*(2^(m+1)+b)',product_residue_class_bases=len(small_bases),
  exact_9_by_9_product_determinants=products,
  scope='Rank 3 of homogenized actual g-index vectors in every dyadic tail; product independence forces bilinear coefficient equality. No coupled-subset or occupancy exclusion.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
