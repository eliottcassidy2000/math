"""Independent integer-only verification; does not import SymPy."""
import argparse,json,math
ap=argparse.ArgumentParser();ap.add_argument('bundle');ap.add_argument('--negative-control',choices=['edge','prime']);args=ap.parse_args()
d=json.load(open(args.bundle));certs=d['prime_certificates'];checked=set()
def need(ok,why):
 if not ok:raise ValueError(why)
def prime(p):
 if p in checked:return
 c=certs[str(p)]
 if p==2:need(c=={'base':2},'prime base')
 else:
  f={int(q):e for q,e in c['p_minus_1_factors'].items()};need(p>2 and all(2<=q<p and e>=1 for q,e in f.items()),'certificate descent')
  need(math.prod(q**e for q,e in f.items())==p-1,'p-1 factor product')
  for q in f:prime(q)
  a=c['witness'];need(pow(a,p-1,p)==1,'Fermat gate')
  for q in f:need(math.gcd(pow(a,(p-1)//q,p)-1,p)==1,'Lucas order gate')
 checked.add(p)
if args.negative_control=='edge':d['orbits'][0]['rows'][0]['s']+=1
if args.negative_control=='prime':certs['9901']['witness']=1
results=[]
for orbit in d['orbits']:
 n=orbit['seed'];seen=set();squares=[]
 for i,r in enumerate(orbit['rows']):
  need(r['i']==i and r['n']==n and n not in seen,'chain or repetition');seen.add(n)
  f={int(p):e for p,e in r['factors'].items()};need(all(e>=1 for e in f.values()),'exponents')
  for p in f:prime(p)
  need(math.prod(p**e for p,e in f.items())==n,'full factor product')
  sig=math.prod(sum(p**j for j in range(e+1)) for p,e in f.items());s=sig-n
  need(r['sigma']==sig and r['s']==s,'aliquot identity')
  g=math.gcd(n,s);need(r['gcd']==g and r['primitive']==[n//g,s//g,sig//g],'normalization')
  sq=math.isqrt(sig)**2==sig;need(r['sigma_square']==sq,'square indicator')
  if sq and s>0:squares.append([n,s])
  if n<=1000001:
   direct=sum(v for k in range(1,math.isqrt(n)+1) if n%k==0 for v in ({k,n//k}-{n}));need(direct==s,'divisor-pair check')
  if f.get(2)==2 and 7 in f:need(s%7==0 and (s%8==4)==(sig%8==0),'valuation rule')
  n=s
  if orbit['seed']==1001000 and i==1106:need(n==34824938264087318322868834739622267021533597146401245964,'first frontier')
  if orbit['seed']==1001000 and i==1136:need(n==120743803148050878826315520222770879154553267874534033029280572,'second frontier')
 need(n==orbit['last_value'],'last value')
 if orbit['status']=='certified_prefix_only':need(n not in seen,'endpoint novelty')
 expected={1000001:(11,0),1001:(5,0),1001000:(1144,3102001427119207440325481174042778229464844199474079305061910724)}[orbit['seed']]
 need((len(orbit['rows']),n)==expected,'frozen count and endpoint')
 need(orbit['status']==('terminal' if n==0 else 'certified_prefix_only'),'status')
 need(squares==([[343,57]] if orbit['seed']==1001 else []),'square-sum census')
 if orbit['seed']==1001000:
  rows=orbit['rows']
  need(all(r['n']%2==0 for r in rows),'observed even prefix')
  for lo,hi in ((709,824),(852,1143)):
   need(all(r['factors'].get('2')==2 and int(r['factors'].get('7',0))>=1 for r in rows[lo:hi+1]),'valuation run')
  need(rows[825]['factors'].get('2')==4,'valuation escape')
 results.append({'seed':orbit['seed'],'edges':len(orbit['rows']),'last_value':n,'square_sum_edges':squares})
need({r['seed'] for r in results}=={1000001,1001,1001000} and len(results)==3,'requested seeds')
R=10071433204932491689368445370268760485275468180110647094356853
need(4*7*11*R==3102001427119207440325481174042778229464844199474079305061910724,'final partial factorization')
need(pow(2,R-1,R)==1059352591669555434863290548381996608006601416003032122563405!=1,'final cofactor composite')
need(5511013429899318414019*15231426307036603838278075031==83940594934600503096994848435730837699779204259589,'first ECM split')
need(7612737786472819497221147*14132068366536545816669370486649==107583730854949977748239831120755574285637006670563966403,'second ECM split')
print(json.dumps({'results':results,'certified_primes':len(checked)},indent=2))
