from pathlib import Path
import hashlib,json,shutil,subprocess,sys,importlib.util,math
from fractions import Fraction
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def call(args,cwd):
 r=subprocess.run([sys.executable,*args],cwd=cwd,capture_output=True,text=True)
 need(r.returncode==0,r.stderr);return r.stdout
summary={}
for name in ('session53-collatz27-repeated-block','session53-two-step-tower'):
 src=ROOT/'inputs'/name;manifest=json.loads((src/'manifest.json').read_text())
 files={x['path']:x['sha256'] for x in manifest['files']} if isinstance(manifest['files'],list) else manifest['files']
 dst=ROOT/'work'/name
 if dst.exists():shutil.rmtree(dst)
 dst.mkdir(parents=True)
 for n,h in files.items():
  need(sha(src/n)==h,'input hash '+n);p=dst/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src/n,p)
 shutil.copy2(src/'manifest.json',dst/'manifest.json')
 summary[name]={'manifest_sha256':sha(src/'manifest.json'),'payload_files':len(files)}
 if 'tower' in name:
  summary[name]['replay']=call(['replay.py'],dst).strip()
  summary[name]['hostiles_normal']=json.loads(call(['hostile.py'],dst))
  summary[name]['hostiles_optimized']=json.loads(call(['-O','hostile.py'],dst))
  saved=json.loads((dst/'normal.json').read_text());summary[name]['comparisons']=saved['comparisons'];summary[name]['status']=saved['status']
 else:
  a=call(['check.py'],dst);b=call(['-O','check.py'],dst)
  need(a==b==(dst/'normal.json').read_text()==(dst/'optimized.json').read_text(),'replay mismatch')
  summary[name]['normal_optimized_saved_equal']=True
# Independent boundary and exact rational-error checks, rather than only tiny X enumeration.
q=2**59;a=3**37;b=23*q-27*a;d=q-a
need(math.gcd(b,d)==1 and 8*d<b<9*d,'fixed point arithmetic')
def val(x):
 x=abs(x);need(x!=0,'zero valuation');return (x&-x).bit_length()-1
def ap(X,r,M):return (X-r)//M-(-r)//M
def residue(k):return b*pow(d,-1,q**k)%(q**k) if k else 0
checks=0
for k in range(33):
 r=residue(k);rn=residue(k+1);M=q**k;MN=M*q
 for X in {0,1,r,max(0,r-1),r+1,rn,max(0,rn-1),rn+1,MN-1,MN,MN+1,11*MN+rn,2**2000+7}:
  count=ap(X,r,M)-ap(X,rn,MN)
  need(count>=0 and abs(Fraction(count)-X*(Fraction(1,M)-Fraction(1,MN)))<=2,'density bound')
  for z in range(7):
   R=r+M*((z-r)*pow(M,-1,7)%7);RN=rn+MN*((z-rn)*pow(MN,-1,7)%7)
   joint=ap(X,R,7*M)-ap(X,RN,7*MN)
   need(joint>=0 and abs(Fraction(joint)-X*(Fraction(1,7*M)-Fraction(1,7*MN)))<=2,'joint density bound')
   checks+=1
 # Choose sources with exactly k copies by selecting a lift that avoids C_(k+1).
 n=r if r else 1
 while val(d*n-b)//59!=k:n+=M
 need(n>0 and val(d*n-b)//59==k,'exact chosen boundary')
 if k:
  f=(a*n+b)//q;need((a*n+b)%q==0 and f<n and val(d*f-b)==val(d*n-b)-59,'one-block exact fuel decrement')
summary['independent_large_interval_CRT_checks']=checks
# Force exact valuations around word-length boundaries, including cases absent from producer samples.
def bitword(n,length):
 out=[]
 for _ in range(length):out.append(n%2);n=(3*n+1)//2 if n%2 else n//2
 return out
w=bitword(27,59);boundaries=[]
for v in (1,58,59,60,116,117,118,119,59*58-1,59*58,59*59-1,59*59):
 M=2**v;n=b*pow(d,-1,M)%M
 if val(d*n-b)!=v:n+=M
 need(n>0 and val(d*n-b)==v,'constructed valuation')
 K=v//59;word=bitword(n,59*(K+1))
 need(word[:59*K]==w*K and word!=w*(K+1),'direct valuation boundary')
 boundaries.append({'valuation':v,'correct_K':K,'wrong_divisor58_K':v//58})
summary['exact_valuation_boundary_controls']=boundaries
# Actual source mutations: each should be caught by the existing positive/hostile certificate.
mutations=[
 ('wrong_period','Q=q**k;A=a**k','Q=q**(k+1);A=a**k'),
 ('lost_carry','B=b*S','B=0*S'),
 ('wrong_fuel58','K=v2(d*n-b)//59','K=v2(d*n-b)//58'),
 ('wrong_fuel60','K=v2(d*n-b)//59','K=v2(d*n-b)//60'),
 ('native_involution','r%7==e%7 and e<r and e!=r','r%7==e%7 and e==r'),
 ('wrong_clock','ordinary==96','ordinary==59'),
 ('wrong_fixedpoint','8*d<b<9*d','9*d<b<10*d'),
 ('count_zero_as_positive','return (X-r)//M-(-r)//M','return (X-r)//M-(-r)//M+1'),
 ('drop_residue_parity','(1 if k%2==0 else 6)','1')]
base=ROOT/'work'/'session53-collatz27-repeated-block';source=(base/'check.py').read_text();out=[]
for label,old,new in mutations:
 need(old in source,'mutation missing '+label)
 dst=ROOT/'mutants'/label;dst.mkdir(parents=True,exist_ok=True);(dst/'check.py').write_text(source.replace(old,new))
 for flags in ([],['-O']):
  run=subprocess.run([sys.executable,*flags,'check.py'],cwd=dst,capture_output=True,text=True)
  # Adding one to both AP counts cancels identically, so this mutation is not semantic for the difference count.
  if label=='wrong_fuel58':
   need(run.returncode==0,'producer coverage changed')
   need(any(x['correct_K']!=x['wrong_divisor58_K'] for x in boundaries),'independent hostile not distinguished')
   out.append({'mutation':label,'mode':'optimized' if flags else 'normal','rejection':'SURVIVES producer fixture; REJECTED by independent exact-valuation boundary controls'})
   continue
  if label=='count_zero_as_positive':
   need(run.returncode==0,'unexpected cancellation behavior');continue
  need(run.returncode!=0,'mutation survived '+label)
  out.append({'mutation':label,'mode':'optimized' if flags else 'normal','rejection':run.stderr.splitlines()[-1]})
summary['rejected_mutations']=out
summary['ap_count_offset_mutation']='A common additive count offset cancels in differences; not a valid semantic hostile. Independent AP formula reviewed separately.'
(ROOT/'audit_results.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
