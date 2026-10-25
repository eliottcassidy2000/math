import pathlib,hashlib,json,subprocess,tempfile,sys,math,itertools
ROOT=pathlib.Path(__file__).resolve().parent
P=ROOT/'inputs/session54-paley-endpoint-observer'
def need(c,m):
 if not c:raise RuntimeError(m)
def validate(p):
 pins=json.loads((p/'MANIFEST.json').read_text());wanted={x['path'] for x in pins}|{'MANIFEST.json'}
 actual={str(x.relative_to(p)) for x in p.rglob('*') if x.is_file()}
 need(actual==wanted,'closed file inventory')
 for x in pins:need(hashlib.sha256((p/x['path']).read_bytes()).hexdigest()==x['sha256'],'manifest pin '+x['path'])
 for x in json.loads((p/'sources.json').read_text()):need(hashlib.sha256((p/x['copy']).read_bytes()).hexdigest()==x['sha256'],'source pin')
validate(P)
with tempfile.TemporaryDirectory() as td:
 binary=pathlib.Path(td)/'independent'
 subprocess.run(['g++','-std=c++17','-O2',str(ROOT/'independent_count.cpp'),'-o',str(binary)],check=True,capture_output=True)
 rows=subprocess.run([str(binary)],check=True,capture_output=True,text=True,timeout=120).stdout
 need(rows==(P/'counts.txt').read_text(),'independent endpoint DP mismatch')
 for opt in [[],['-O']]:
  replay=subprocess.run([sys.executable,*opt,str(P/'verify.py')],check=True,capture_output=True,text=True,timeout=120)
  need(replay.stdout==(P/'receipt.json').read_text(),'producer normal/-O receipt')
records=[]
for line in rows.splitlines():
 q,*v=map(int,line.split());r=(q-1)//2;R={k*k%q for k in range(1,q)};a=v[1];b=v[next(x for x in range(1,q) if x not in R)];s=a+b;delta=a-b;t=r*s;N=(s*s+q*delta*delta)//4;M=t*t+s*t+N
 need(v==[0]+[a if x in R else b for x in range(1,q)],'endpoint channel')
 need(s%2==delta%2==1 and delta!=0,'parity')
 need(M==q*(q*s*s+delta*delta)//4,'resultant modulus')
 need(b%q==1,'translation congruence')
 need(math.factorial(q)<2**64,'count arithmetic bound')
 # Basis 1,E,E^2 in Z x Z[lambda] has columns (1,1,0),(t,0,1),(t^2,-N,-s).
 determinant=N+s*t+t*t
 need(determinant==M and M>1,'integral product hostile')
 # Arbitrary polynomial is evaluated with lambda^2=-s lambda-N.
 for coeffs in itertools.product(range(-1,2),repeat=5):
  u=sum(c*t**i for i,c in enumerate(coeffs));A=B=0
  for c in reversed(coeffs):A,B=c-N*B,A-s*B
  need((u-A-B*t)%M==0,'fiber congruence')
 # Conversely A+BX+kQ(X) attains all triples with difference kM.
 for A,B,k in itertools.product(range(-2,3),repeat=3):
  u=A+B*t+k*M
  need(u-(A+B*t)==k*M,'fiber sufficiency')
 need(s*s+s*s+N>0,'scalar does not solve quadratic')
 need(-q*delta*delta<0,'irreducible')
 if q==19:need(a%q==3 and delta%q==2,'q19 extrapolation hostile')
 if q==3:need(a==0 and b==1,'strict positivity hostile')
 records.append(dict(q=q,a=a,b=b,H=q*r*s,C=r*b,s=s,delta=delta,N=N,M=M,conductor=abs(delta),discriminant=-q*delta*delta,independent_dp_bytes=(1<<q)*q*8))
print(json.dumps(dict(status='PASS',producer_manifest_sha256=hashlib.sha256((P/'MANIFEST.json').read_bytes()).hexdigest(),records=records,scope='Prime q theorem; independent exact DP q=3,7,11,19; no modular coefficient law; integral order not direct product'),indent=2))
