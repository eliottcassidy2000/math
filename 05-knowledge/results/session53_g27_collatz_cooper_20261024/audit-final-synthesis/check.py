from pathlib import Path
from math import comb
import hashlib,json,re
P=Path(__file__).resolve().parent; sha=lambda b:hashlib.sha256(b).hexdigest()
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def entries(p):
 d=json.loads(p.read_text());d=d.get('files',d)
 return {x['path']:x['sha256'] for x in d} if isinstance(d,list) else d
pins=json.loads((P/'SOURCE_PINS.json').read_text())
for name,pin in pins.items():
 d=P/'sources'/name;m=d/pin['manifest'];need(sha(m.read_bytes())==pin['sha256'],'pin '+name);e=entries(m)
 need({x.relative_to(d).as_posix() for x in d.rglob('*') if x.is_file()}==set(e)|{m.name},'source closure '+name)
 for f,h in e.items():need(sha((d/f).read_bytes())==h,'source hash '+f)
link_count=0
for name in ['README.md','SYNTHESIS.md']:
 text=(P/('REVIEWED_'+name)).read_text()
 for link in re.findall(r'\]\(([^)]+)\)',text):
  target=P/('REVIEWED_'+link) if link in ['README.md','SYNTHESIS.md'] else P/'sources'/link
  need(target.is_file(),'link '+link);link_count+=1
 need('v 2(' not in text,'ambiguous valuation typography')
need('e94dfc15927a7548d8e942cdf74aa94b030f660e' in (P/'REVIEWED_README.md').read_text(),'commit identity')
c=[];f=[];weights=[1,2,4,3,3,6,0,0,0]
for n in range(101):
 f.append(sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range(n+1)))
 if not n:c.append(1)
 else:
  v=f[n]-sum(c[k]*c[n-k] for k in range(1,n));need(v%2==0,'root integrality');c.append(v//2)
 need(c[n]%4==(1 if n==0 else 2 if n&(n-1)==0 else 0),'dyadic')
 v=n;z=1
 while v:z=z*weights[v%9]%9;v//=9
 need(c[n]%9==z,'base9')
need(f[3]-f[1]==756 and 756%81==27,'excellent sharpness')
need(f[5]-48*f[1]==273312 and 273312%27==18,'naive coordinate')
need(c[9]==26337699740 and (c[9]-c[1])%27==18,'two step sharpness')
primes=[p for p in range(2,201) if all(p%d for d in range(2,int(p**.5)+1))]
holdouts=sum(40000//(p*p)-1 for p in primes if p>2 and p*p*2<=40000)
need(holdouts==7995,'comparison count');need(sum(p>2 and p*p*2<=40000 for p in primes)==33,'prime count')
tower=0
for p in primes:
 power=p*p
 while power<=20000:tower+=20000//power;power*=p
need(tower==15407,'tower comparisons')
q=2**59;a=3**37;b=23*q-27*a;d=q-a
need(8*d<b<9*d and d%2==1 and a<q,'fixed point')
x=27;word=[]
for j in range(59):
 word.append(x%2);x=(3*x+1)//2 if x%2 else x//2
need(x==23 and sum(word)==37,'clock')
need(a*pow(q,-1,7)%7==6 and b*pow(q,-1,7)%7==1,'residue action')
boundaries=[0,1,57,58,59,60,116,117,118,119,176,177,235,236]
for v in boundaries:
 mod=2**(v+1);n=((b+2**v)*pow(d,-1,mod))%mod
 if n==0:n+=mod
 delta=abs(d*n-b);valuation=(delta&-delta).bit_length()-1
 need(valuation==v,'exact valuation construction');x=n;k=0
 for block in range(v//59+1):
  bits=[]
  for j in range(59):bits.append(x%2);x=(3*x+1)//2 if x%2 else x//2
  if bits!=word:break
  k+=1
 need(k==v//59,'repeat fuel')
# Independent union-bound identity on finite coefficient range; density proof is in REPORT.
for n in range(101):need((c[n]%36==0)==(c[n]%4==0 and c[n]%9==0),'coprime intersection')
# Shared determinant/orientation sign, not an operator identification.
chi=lambda x:0 if x%7==0 else 1 if x%7 in (1,2,4) else -1
for residue in range(1,7):
 determinant=1 if residue in (1,2,4) else -1
 need(determinant==chi(residue),'Cartier determinant sign')
 for x in range(7):
  for y in range(7):
   if x!=y:need(chi(residue*(y-x))==determinant*chi(y-x),'Paley difference orientation')
need(len({7*x%7 for x in range(7)})==1,'ramified nonpermutation')
print(json.dumps({'status':'PASS','reviewed_hashes':{x:sha((P/('REVIEWED_'+x)).read_bytes()) for x in ['README.md','SYNTHESIS.md']},'source_packets':len(pins),'relative_link_uses':link_count,'independent_binomial_root_through':100,'next_digit_comparisons_outside_calibration':holdouts,'next_digit_odd_primes':33,'tower_comparisons_counted':tower,'exact_fuel_valuation_boundaries':boundaries,'sharpness':{'excellent_mod81':756,'naive_mod27':273312}},indent=2))
