from pathlib import Path
from math import comb
import hashlib,json,re
P=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def entries(p):
 d=json.loads(p.read_text());d=d.get('files',d)
 return {x['path']:x['sha256'] for x in d} if isinstance(d,list) else d
pins=json.loads((P/'SOURCE_PINS.json').read_text())
for name,pin in pins.items():
 d=P/'sources'/name;m=d/pin['manifest'];need(sha(m.read_bytes())==pin['sha256'],'pin '+name);e=entries(m)
 need({x.relative_to(d).as_posix() for x in d.rglob('*') if x.is_file()}==set(e)|{m.name},'closed source '+name)
 for f,h in e.items():need(sha((d/f).read_bytes())==h,'payload '+f)
links=re.findall(r'\]\(([^)]+)\)',(P/'SUMMARY.md').read_text())
for link in links:need((P/'sources'/link).is_file(),'summary link '+link)
a=entries(P/'lineage/first/MANIFEST.json');b=entries(P/'lineage/hardened/MANIFEST.json');need({k:v for k,v in a.items() if k!='replay.py'}=={k:v for k,v in b.items() if k!='replay.py'},'lineage mathematical hashes')
for tag,e in [('first',a),('hardened',b)]:need(sha((P/'lineage'/tag/'replay.py').read_bytes())==e['replay.py'],'replay hash')
weights=[1,2,4,3,3,6,0,0,0]; c=[]
for n in range(101):
 f=sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range(n+1))
 if n==0:c.append(1)
 else:
  q=f-sum(c[k]*c[n-k] for k in range(1,n));need(q%2==0,'integral root');c.append(q//2)
 v=n;z=1
 while v:z=z*weights[v%9]%9;v//=9
 need(c[n]%9==z,'digit '+str(n))
need(c[9]==26337699740 and (c[9]-c[1])%27==18,'mod27 obstruction')
counts=[0,1]+[0]*7;h=[2,1,-1,-2,-1,1]
for r in range(13):
 need(counts[0]==9**r-(r+1)*3**r,'zero count')
 need(counts[3]==(r*3**r+r)//2 and counts[6]==(r*3**r-r)//2,'defect counts')
 if r:
  for j in range(6):need(6*counts[pow(2,j,9)]==3**r+(-1)**j+2**r*h[(r-j)%6],'unit counts')
 nxt=[0]*9
 for residue,k in enumerate(counts):
  for w in weights:nxt[residue*w%9]+=k
 counts=nxt
print(json.dumps({'status':'PASS','summary_sha256':sha((P/'SUMMARY.md').read_bytes()),'closed_source_packets':len(pins),'summary_relative_link_uses':len(links),'coefficient_binomial_convolution_range':[0,100],'distribution_digit_places':[0,12],'lineage_changed_payloads':['replay.py'],'source_pins':pins},indent=2))
