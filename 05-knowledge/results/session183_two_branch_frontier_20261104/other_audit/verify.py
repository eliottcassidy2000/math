import hashlib,json,itertools
from pathlib import Path
OUT=Path(__file__).parent
TARGET=OUT.parent/'other'
def need(x,s):
 if not x:raise RuntimeError(s)
def j_by_covers(f):
 m=f.bit_length();q=(1<<m)-1;z=0
 for left in range(1,q):
  right=q^left
  z+=bool((f>>(left.bit_length()-1)&1) and (f>>((right&-right).bit_length()-1)&1))
 return z
# Literal cover evaluator independent of both target's low-bit recursion and its high-bit scan.
j=[j_by_covers(f) for f in range(3*512+3)]
def D(c,n):return 3*j[n]-j[3*n+c]
for n in range(256):
 p=n%2
 need([D(0,2*n),D(0,2*n+1),D(1,2*n),D(1,2*n+1),D(2,2*n),D(2,2*n+1)]==[D(0,n),D(1,n)+4*p-2,D(0,n)-3*n-p,D(2,n)+3*n+3*p,D(1,n),D(2,n)+2*p-2],('six carries',n))
for n in range(512):
 need(-2*n<=D(0,n)<=2*n and -3*n<=D(1,n)<=4*n and -2*n<=D(2,n)<=2*n,('bounds',n))
 if n:need(D(0,n)!=4 and (D(0,n)==6)==(n//(n&-n)==7),('classification',n))
need([D(0,4*l+3) for l in range(5)]==[2,6,8,12,8],'five base cases')
normal=(OUT/'normal.json').read_bytes();optimized=(OUT/'optimized.json').read_bytes()
need(normal==optimized,'normal and optimized match')
need(json.loads(normal)==json.loads((TARGET/'verification.json').read_text()),'published receipt matches rerun')
result={'status':'PASS','literal_cover_masks':len(j),'independent_six_carry_checks':256,'independent_bounds_checks':512,'independent_classification_checks':511,'rerun_normal_optimized_identical':True,'published_receipt_matches':True,'reviewed_files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [TARGET/'PROOF.md',TARGET/'verify.py',TARGET/'verification.json']}}
(OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
