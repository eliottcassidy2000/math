from core import *
from pruning import prune
from pathlib import Path
import itertools,json,hashlib,random,subprocess,tempfile,os
P=Path(__file__).resolve().parent

def check_pruning():
 controls=[];digest=hashlib.sha256();total=0
 for k,depth in [(0,5),(1,4),(2,3),(3,2)]:
  n=max(k,1);seed=[sum(1<<j for j in range(i+1,n))for i in range(n)];ports=list(range(k))
  alphabet=list(itertools.product(('R','L'),range(1<<k),range(k+1)));count=0
  for word in itertools.product(alphabet,repeat=depth):
   for mask in range(1<<depth):
    keep=set(range(n))|{n+j for j in range(depth)if mask>>j&1}
    rows,pp,w,shadow=prune(seed,ports,word,keep)
    require(len(w)==mask.bit_count(),'exact number of kept births')
    digest.update(str((rows,pp,w,sorted(shadow.items()))).encode());count+=1
  controls.append({'ports':k,'depth':depth,'induced_subsequence_tests':count});total+=count
 rng=random.Random(197541);long=0
 for j in range(100):
  word=[(rng.choice(('R','L')),rng.randrange(8),rng.randrange(4))for _ in range(12)]
  for _ in range(30):
   mask=rng.randrange(1<<12);keep=set(range(8))|{8+i for i in range(12)if mask>>i&1}
   rows,pp,w,shadow=prune(SEED,PORTS,word,keep)
   require(len(w)==mask.bit_count(),'long pruning length');long+=1
   digest.update(str((rows,pp,w,sorted(shadow.items()))).encode())
 return {'exhaustive':controls,'exhaustive_total':total,'long_controls':long,'data_sha256':digest.hexdigest()}

def bound_controls():
 rng=random.Random(197179);records=[]
 def omega(n):
  ans=0;d=2
  while d*d<=n:
   while n%d==0:n//=d;ans+=1
   d+=1
  return ans+(n>1)
 def blog(n):
  t=0;p=1
  while p*5<=n**3:p*=5;t+=1
  return t
 for j in range(200):
  seed=[2,4,1];ports=[0,1,2];rows=seed;pp=ports;word=[]
  for _ in range(5):
   op=(rng.choice(('R','L')),rng.randrange(8),rng.randrange(4));word.append(op)
   rows,pp=birth(rows,pp,*op)
  reach=[r|1<<i for i,r in enumerate(rows)]
  for k in range(len(rows)):
   for i in range(len(rows)):
    if reach[i]>>k&1:reach[i]|=reach[k]
  keep=set(range(3));nontrivial=set()
  for i in range(len(rows)):
   C={v for v in range(len(rows))if (reach[i]>>v&1)and(reach[v]>>i&1)}
   if len(C)>1:nontrivial|=C
  keep|=nontrivial;N=hp(rows);B=omega(N)+blog(N)
  out,qp,ww,shadow=prune(seed,ports,word,keep)
  require(hp(out)==N,'SCC pruning preserves H')
  require(len(nontrivial)<=B and len(ww)<=B-3,'strong-seed finite depth bound')
  records.append([N,len(nontrivial),len(ww),B])
 return {'cases':len(records),'data_sha256':hashlib.sha256(str(records).encode()).hexdigest()}

def python_census2():
 frontier=[(SEED,PORTS,179)];sizes=[1]
 for depth in range(2):
  nxt=[]
  for rows,ports,H in frontier:
   for side in ('R','L'):
    for mask in range(8):
     rr,_=birth(rows,ports,side,mask,3);h=hp(rr)
     require(h!=541,'independent short scalar target absence')
     if h>541:continue
     for drop in range(4):
      qq=ports+[len(rows)];del qq[drop];nxt.append((rr,qq,h))
  frontier=nxt;sizes.append(len(frontier))
 require(sizes==[1,36,440],'independent short census sizes')
 return sizes

def main():
 require(all(541%d for d in range(2,24)),'target prime')
 require(5**12>541**3,'strong target component order at most12')
 require(len(SEED)==8 and strong(SEED) and hp(SEED)==179,'strong seed')
 plus2,pp=birth(SEED,PORTS,'R',7,3)
 require(hp(plus2)==539 and pp==PORTS,'positive physical next clock step')
 pc=check_pruning();small=python_census2();bc=bound_controls()
 with tempfile.TemporaryDirectory(prefix='prune-census-',dir=P)as td:
  exe=Path(td)/'census';bank=Path(td)/'bank.jsonl'
  flags=['-O0','-D_GLIBCXX_ASSERTIONS']if os.environ.get('PRUNING_CPP_DEBUG')else['-O3']
  subprocess.run(['g++','-std=c++17',*flags,str(P/'census.cpp'),'-o',str(exe)],check=True)
  r=subprocess.run([str(exe),str(bank)],check=True,text=True,capture_output=True);c=json.loads(r.stdout)
  require(c['accepted_word_prefixes']==[1,36,440,5716,79768],'full exact census size')
  require(all('541'not in h for h in c['histograms']),'complete target absence')
  c['bank_sha256']=hashlib.sha256(bank.read_bytes()).hexdigest()
  require(sum(c['accepted_word_prefixes'])==sum(1 for _ in bank.open()),'complete bank serialization')
 return {'status':'PASS','seed':SEED,'ports':PORTS,'target':541,'prime':True,
         'maximum_strong_order':12,'maximum_kept_births':4,'pruning':pc,
         'independent_python_census2':small,'bound_controls':bc,'census':c}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
