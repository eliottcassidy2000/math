from core import *
from pathlib import Path
from collections import Counter
import itertools,json,hashlib,subprocess,tempfile,os
P=Path(__file__).resolve().parent
REBASE=[2,64,11,51,39,135,61,95]
RP=[0,1,2,5]
def canonical(rows):
 n=len(rows);best=None
 for tail in itertools.permutations(range(8,n)):
  p=tuple(range(8))+tail
  out=tuple(sum(((rows[u]>>v)&1)<<j for j,v in enumerate(p))for u in p)
  if best is None or out<best:best=out
 return best

def literal(rows):
 return sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:]))for p in itertools.permutations(range(len(rows))))

def main():
 require(hp(SEED)==179 and strong(SEED),'original strong seed')
 require(all(541%d for d in range(2,24)),'prime541')
 require(5**12>541**3>=5**11,'exact strong order cap')
 with tempfile.TemporaryDirectory(prefix='supergraph-census-',dir=P)as td:
  exe=Path(td)/'census';bf=Path(td)/'bank.jsonl';flags=['-O0','-D_GLIBCXX_ASSERTIONS']if os.environ.get('SUPERGRAPH_CPP_DEBUG')else['-O3']
  subprocess.run(['g++','-std=c++17',*flags,str(P/'census.cpp'),'-o',str(exe)],check=True)
  run=subprocess.run([str(exe),str(bf)],check=True,text=True,capture_output=True);stats=json.loads(run.stdout)
  bank=[json.loads(line)for line in bf.read_text().splitlines()]
  require(dict(Counter(z['n']for z in bank))=={8:1,9:33,10:68,11:105,12:144},'complete bank sizes')
  seen=set()
  for z in bank:
   n=z['n'];r=z['rows'];H=z['H'];key=(n,tuple(r));require(key not in seen,'no duplicate');seen.add(key)
   require(len(r)==n and all(0<=x<1<<n and not x>>i&1 for i,x in enumerate(r)),'valid rows')
   require(all(((r[i]>>j)&1)+((r[j]>>i)&1)==1 for i in range(n)for j in range(i+1,n)),'tournament')
   require([x&255 for x in r[:8]]==SEED,'fixed induced seed')
   require(canonical(r)==tuple(r),'exact suffix canonicalization')
   require(hp(r)==H<=541 and H!=541,'independent graph count')
  control=set()
  for mask in range(256):
   r=SEED[:]+[mask]
   for j in range(8):
    if not mask>>j&1:r[j]|=256
   if hp(r)<=541:control.add(tuple(r))
  require(control=={tuple(z['rows'])for z in bank if z['n']==9},'independent complete one-vertex extension bank')
  stats['bank_sha256']=hashlib.sha256(bf.read_bytes()).hexdigest()
  stats['independent_bank_checks']=len(bank)
 plus2,p2=birth(REBASE,RP,'L',4,4);plus4,p4=birth(REBASE,RP,'R',11,4)
 require(p2==p4==RP,'same retained four-port interface')
 witnesses=[]
 for r,h in [(REBASE,179),(plus2,539),(plus4,541)]:
  require(strong(r)and hp(r)==h and literal(r)==h,'independent physical rebasing witnesses')
  witnesses.append({'rows':r,'H':h,'strong':True})
 return {'status':'PASS','seed':SEED,'target':541,'maximum_strong_order':12,
         'census':stats,'rebased_ports':RP,'rebase_witnesses':witnesses,
         'rebase_words':[['L',4,4],['R',11,4]]}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
