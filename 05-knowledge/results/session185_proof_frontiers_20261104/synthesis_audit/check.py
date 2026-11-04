import re,json,hashlib
from pathlib import Path
from fractions import Fraction
P=Path(__file__).parent;B=P.parent;note=B/'synthesis'/'NOTE.md'
def need(b,m):
 if not b:raise RuntimeError(m)
s=note.read_text();links=[]
for link in re.findall(r'\]\(([^)]+)\)',s):
 if '://' not in link:
  f=(note.parent/link).resolve();need(f.is_file(),'missing sibling '+link);links.append(link)
w={x['address']:x for x in json.loads((B/'closure'/'witnesses.json').read_text())}
sig={'2':(21,-13),'4':(21,-11)}
for a,v in w.items():sig[a]=(v['sig'][2],v['sig'][0]-v['sig'][4])
edges=[('2','24',4,-2),('4','44',4,-3),('24','242',2,0),('24','244',4,-1),('44','442',2,1),('44','444',4,2)]
for a,b,d,q in edges:
 f,c=sig[a];F,C=sig[b];need(Fraction(2)**q*F==3*f and C==3*c+d,'edge '+b)
clocks={'2':4,'4':4,'24':2,'44':1,'242':2,'244':1,'442':2,'444':3}
counts={a:c+(1<<clocks[a])*f for a,(f,c) in sig.items()}
need(counts=={'2':323,'4':325,'24':973,'44':979,'242':2921,'244':2923,'442':2939,'444':2941},'count paths')
for a in ('24','44'):
 z=w[a]['sig'];c=z[0]-z[4];need(c<=-3 and min(z[2],z[3])>-2*c-2,'context hypothesis')
for r,expected in [(1,[27,29,31]),(2,[69,71,73]),(3,[177,179,181])]:need([4*3**r+9*2**r-3+2*j for j in range(3)]==expected,'gate interval')
for r,expected in ((4,469),(6,3493)):
 H=4*3**r+9*2**r+1;need(H==expected and H%7==0 and H%49!=0,'gate valuation')
need(pow(2,42,49)==pow(3,42,49)==1,'gate period')
for e in range(2,12):need(Fraction(1,3*7**(e-2))-Fraction(1,3*7**(e-1))==Fraction(2,7**(e-1)),'density')
need(7*(27*2**35-13)==6493990551461,'large prime-family count')
files=[note]+[B/d/'PROOF.md' for d in ('closure','mining','carriers','alternative','prime_lanes')]
files +=[B/(d+'_audit')/'AUDIT.md' for d in ('closure','mining','carriers','alternative','prime_lanes')]
out={'status':'PASS','relative_links_checked':links,'count_paths':counts,'audited_sha256':{str(f.relative_to(B)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},'scope':'Synthesis cross-packet consistency and arithmetic; remote pinned links not refetched.'}
(P/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
