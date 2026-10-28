"""Regenerate exact aliquot certificates with installed SymPy 1.14.0.
Usage: python3 generate_aliquot_orbits.py --output /tmp/aliquot.json
The four exact factor hints avoid repeating the two explicit ECM discovery searches.
Ordinary SymPy factorint may still use ECM internally.
"""
import argparse,json,math,signal
import sympy as sp
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--seconds',type=int,default=60);args=ap.parse_args()
if args.seconds<=0:ap.error('--seconds must be positive')
if not hasattr(signal,'SIGALRM'):ap.error('POSIX SIGALRM support is required')
HINTS=(5511013429899318414019,15231426307036603838278075031,7612737786472819497221147,14132068366536545816669370486649)
CHECKPOINTS={1107:34824938264087318322868834739622267021533597146401245964,1137:120743803148050878826315520222770879154553267874534033029280572,1144:3102001427119207440325481174042778229464844199474079305061910724}
certs={}
def fail_timeout(*_):raise TimeoutError('single factorization budget exceeded; no uncertified edge emitted')
signal.signal(signal.SIGALRM,fail_timeout)
def factor(n):
 original=n;f={}
 for p in HINTS:
  while n%p==0:f[p]=f.get(p,0)+1;n//=p
 signal.alarm(args.seconds)
 try:
  for p,e in sp.factorint(n).items():f[int(p)]=f.get(int(p),0)+int(e)
 finally:signal.alarm(0)
 if math.prod(p**e for p,e in f.items())!=original:raise ValueError('factor product')
 for p in f:prime(p)
 return f
def prime(p):
 if str(p) in certs:return
 if p==2:certs[str(p)]={'base':2};return
 if not sp.isprime(p):raise ValueError('composite factor')
 f=factor(p-1)
 for a in range(2,10000):
  if pow(a,p-1,p)==1 and all(math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in f):
   certs[str(p)]={'p_minus_1_factors':f,'witness':a};return
 raise ValueError('no Lucas witness found')
orbits=[]
for seed,limit in ((1000001,11),(1001,5),(1001000,1144)):
 n=seed;rows=[];seen=set()
 while n and len(rows)<limit:
  if n in seen:raise ValueError('unexpected cycle')
  seen.add(n);f=factor(n);sigma=math.prod(sum(p**j for j in range(e+1)) for p,e in f.items());s=sigma-n;g=math.gcd(n,s)
  rows.append({'i':len(rows),'n':n,'factors':f,'sigma':sigma,'s':s,'gcd':g,'primitive':[n//g,s//g,sigma//g],'sigma_square':math.isqrt(sigma)**2==sigma});n=s
  if seed==1001000 and len(rows) in CHECKPOINTS and n!=CHECKPOINTS[len(rows)]:raise ValueError('checkpoint mismatch')
 if len(rows)!=limit or (seed!=1001000 and n!=0):raise ValueError('length or terminal mismatch')
 orbits.append({'seed':seed,'rows':rows,'last_value':n,'status':'terminal' if n==0 else 'certified_prefix_only'})
with open(args.output,'w') as f:json.dump({'sympy':sp.__version__,'prime_certificates':certs,'orbits':orbits},f)
print('Generated',sum(len(d['rows']) for d in orbits),'certified edges; C stops intentionally before evaluating n_1144.')
