import json,math,itertools,hashlib,subprocess
from pathlib import Path

def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def square(n):return n>=0 and math.isqrt(n)**2==n
base=Path(__file__).resolve().parent.parent
p=next((base/name for name in ('session42-gaussian','gaussian') if (base/name/'MANIFEST.json').exists()),None)
check(p is not None,'missing sibling Gaussian packet')
repo=next((r for r in [base/'session35-publication/clean-repo',*Path(__file__).resolve().parents] if (r/'.git').exists()),None)
check(repo is not None,'missing repository for inherited source pins')
manifest=json.loads((p/'MANIFEST.json').read_text())
for f in manifest['files']:
 b=(p/f['path']).read_bytes();check(len(b)==f['bytes'] and hashlib.sha256(b).hexdigest()==f['sha256'],'file pin')
provenance=json.loads((p/'sources.json').read_text())
for f in provenance['sources']:
 b=subprocess.check_output(['git','-C',str(repo),'show',provenance['commit']+':'+f['path']])
 check(hashlib.sha256(b).hexdigest()==f['sha256'],'source pin')
# Independent Burnside sum, not representative orbit-size classification.
for m in range(251):
 pts=[(x,y) for x in range(-math.isqrt(m),math.isqrt(m)+1) for y in range(-math.isqrt(m),math.isqrt(m)+1) if x*x+y*y<=m]
 fixed=0
 for swap in (0,1):
  for s,t in itertools.product((-1,1),repeat=2):
   fixed+=sum((s*(y if swap else x),t*(x if swap else y))==(x,y) for x,y in pts)
 reps={(max(abs(x),abs(y)),min(abs(x),abs(y))) for x,y in pts}
 check(fixed==8*len(reps),'Burnside independent disk')
# Both directions of quartic realizability, including unrealizable candidate data.
quartic=0
for u in range(101):
 actual={x*x*y*y for x in range(math.isqrt(u)+1) for y in range(x+1) if x*x+y*y==u}
 for v in range(u*u//4+1):
  d=u*u-4*v
  predicted=False
  if square(d):
   r=math.isqrt(d)
   predicted=(u+r)%2==0 and square((u+r)//2) and square((u-r)//2)
  check(predicted==(v in actual),'quartic realizability');quartic+=1
# Finite exact angle harmonic and joint-frame counterexample.
for x,y in itertools.product(range(-20,21),repeat=2):
 u=x*x+y*y;v=x*x*y*y
 check(x**4-6*x*x*y*y+y**4==u*u-8*v,'quartic harmonic')
check((1+1)**2 != (1-1)**2,'addition hostile')
# Conjugation and norm for every mod-three golden element.
for a,b in itertools.product(range(3),repeat=2):
 image=((a+2*b)%3,b)
 conj=((a+b)%3,(-b)%3)
 mapped_conj=((conj[0]+2*conj[1])%3,conj[1])
 check(mapped_conj==(image[0],-image[1]%3),'conjugation bridge')
 check((a*a+a*b-b*b)%3==(image[0]**2+image[1]**2)%3,'norm bridge')
# Primitive normalization and Pell-locus boundary.
normalization=0
for n in range(2,61):
 for m in range(1,n):
  if math.gcd(m,n)!=1:continue
  a,b,c=n*n-m*m,2*m*n,n*n+m*m
  check(math.gcd(math.gcd(a,b),c)==(2 if n%2==m%2 else 1),'primitive content')
  normalization+=1
fib=[1,1]
while fib[-1]<=500:fib.append(sum(fib[-2:]))
expected={(fib[i],fib[i+1]) for i in range(len(fib)-1) if fib[i+1]<=500}
actual={(m,n) for m in range(1,501) for n in range(1,501) if abs(n*n-m*n-m*m)==1}
check(actual==expected,'Pell locus control')
check((1*1-1*1,2*1*1,1*1+1*1)==(0,2,2),'degenerate terminal')
print(json.dumps({'status':'PASS','manifest_sha256':hashlib.sha256((p/'MANIFEST.json').read_bytes()).hexdigest(),'manifest_files':len(manifest['files']),'source_pins':len(provenance['sources']),'independent_burnside_disks':251,'quartic_candidates':quartic,'primitive_normalization_pairs':normalization,'pell_bound':500},indent=2))
